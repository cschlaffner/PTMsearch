from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MSSpectrum

from src.diagnostic_ions.utils import get_modification_unimod_format
from src.mzml_processing.utils import (
    check_ms1_spectrum,
    get_ms2_spectrum_mz,
    validate_collision_energy_ms2_spectrum,
    validate_ms1_spectrum,
)


class ScanWindowSplitting:

    def __init__(
        self,
        lower_collision_energy: Union[float, int],
        higher_collision_energy: Union[float, int],
    ) -> None:
        self.lower_collision_energy = lower_collision_energy
        self.higher_collision_energy = higher_collision_energy

    def _list_spectra_by_ms1_and_mz(
        self, ms1_and_ms2_spectra: List[MSSpectrum], num_ms1_spectra: int
    ) -> pd.DataFrame:
        """Creates a DataFrame where every MS2 spectrum is identified by its center
        m/z and corresponding MS1 spectrum ID. This is required to match higher- and
        lower-energy spectra."""
        num_ms2_spectra = len(ms1_and_ms2_spectra) - num_ms1_spectra
        assert check_ms1_spectrum(
            ms1_and_ms2_spectra[0]
        ), "Scan window list should start with an MS1 spectrum."
        current_ms1_spectrum = ms1_and_ms2_spectra[0]

        spectra_df = pd.DataFrame(
            columns=[
                "ms1_spectrum_id",
                "ms2_spectrum_mz",
                "ms2_spectrum_id",
                "ms2_spectrum",
            ],
            index=range(num_ms2_spectra),
        )

        ms2_spectra_count = 0
        for spectrum in ms1_and_ms2_spectra:
            if check_ms1_spectrum(spectrum):
                current_ms1_spectrum = spectrum
            else:
                spectra_df.loc[ms2_spectra_count] = [
                    current_ms1_spectrum.getNativeID(),
                    get_ms2_spectrum_mz(spectrum),
                    spectrum.getNativeID(),
                    spectrum,
                ]
                ms2_spectra_count += 1

        return spectra_df

    def _list_ms1_spectra_by_id(self, ms1_windows: MSExperiment) -> pd.DataFrame:
        """Helper function because normal MSExperiment does not allow retrieving
        spectra by native ID.
        """
        spectra = ms1_windows.getSpectra()
        spectra_df = pd.DataFrame(
            columns=["ms1_spectrum_id", "ms1_spectrum"], index=range(len(spectra))
        )
        for i, spectrum in enumerate(spectra):
            spectra_df.loc[i] = [
                spectrum.getNativeID(),
                spectrum,
            ]
        spectra_df.set_index("ms1_spectrum_id", inplace=True)
        return spectra_df

    def _validate_spectra_identifiers_matching(
        self,
        ms1_windows_df: pd.DataFrame,
        lower_energy_windows_df: pd.DataFrame,
        higher_energy_windows_df: pd.DataFrame,
    ) -> None:
        assert lower_energy_windows_df[["ms1_spectrum_id", "ms2_spectrum_mz"]].equals(
            higher_energy_windows_df[["ms1_spectrum_id", "ms2_spectrum_mz"]]
        ), "MS1 windows and/or precursor m/z do not match between lower- and higher-energy windows."

        # Using numpy instead of panda equal methods to make it about values only
        assert np.array_equal(
            ms1_windows_df.index,
            lower_energy_windows_df["ms1_spectrum_id"].drop_duplicates(),
        ), "MS1 window list does not match the MS1 windows provided with the lower- and higher-energy windows."

    def _validate_spectra_mslevel_and_collision_energy(
        self,
        ms1_windows: MSExperiment,
        ms1_and_lower_energy_windows: MSExperiment,
        ms1_and_higher_energy_windows: MSExperiment,
    ):
        for spectrum in ms1_windows.getSpectra():
            validate_ms1_spectrum(spectrum)

        for spectrum in ms1_and_lower_energy_windows.getSpectra():
            if not check_ms1_spectrum(spectrum):
                validate_collision_energy_ms2_spectrum(
                    spectrum, self.lower_collision_energy
                )

        for spectrum in ms1_and_higher_energy_windows.getSpectra():
            if not check_ms1_spectrum(spectrum):
                validate_collision_energy_ms2_spectrum(
                    spectrum, self.higher_collision_energy
                )

    def _match_higher_energy_spectra_with_detected_ions(
        self, higher_energy_windows_df: pd.DataFrame, detected_ions_df: pd.DataFrame
    ) -> pd.DataFrame:
        """The data about the detected ions contains only the ID of the higher-energy spectrum.
        To match the ions with the lower-energy spectra, more higher-energy window information
        (center m/z and MS1 spectrum ID) is needed.
        """

        higher_energy_windows_df.set_index("ms2_spectrum_id", inplace=True)

        # Here, we do not need any additional information (e.g. the ion type) anymore
        detected_ions_df = detected_ions_df[
            ["spectrum_id", "letter_and_unimod_format_mod"]
        ].drop_duplicates()

        # Group detected ions by mod combinations
        detected_ions_df = (
            detected_ions_df.groupby("spectrum_id")["letter_and_unimod_format_mod"]
            .apply(list)
            .to_frame()
        )

        higher_energy_windows_detected_ions_df = detected_ions_df.join(
            higher_energy_windows_df, on="spectrum_id", how="right"
        )
        higher_energy_windows_detected_ions_df.loc[
            higher_energy_windows_detected_ions_df[
                "letter_and_unimod_format_mod"
            ].isna(),
            "letter_and_unimod_format_mod",
        ] = ["unmodified"]

        return higher_energy_windows_detected_ions_df

    def _match_lower_and_higher_energy_spectra_with_ions(
        self,
        lower_energy_windows_df: pd.DataFrame,
        higher_energy_windows_detected_ions_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Matches the lower-energy spectra with the information about the detected ions
        containing also center m/z and MS1 spectrum ID for each higher-energy window
        in order to perform the matching.
        """

        # re-index to improve performance during join
        lower_energy_windows_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_mz"], inplace=True
        )
        # TODO: reset index beforehand?
        higher_energy_windows_detected_ions_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_mz"], inplace=True
        )

        windows_with_detected_ions_df = higher_energy_windows_detected_ions_df.join(
            lower_energy_windows_df,
            on=["ms1_spectrum_id", "ms2_spectrum_mz"],
            how="inner",  # although that shouldn't matter here
            lsuffix="_higher_energy",
        )
        return windows_with_detected_ions_df

    def _group_windows_by_mods(
        self,
        ms1_windows_df: pd.DataFrame,
        windows_with_detected_ions_df: pd.DataFrame,
        mod_combinations_to_search: List[Tuple[str]],
    ) -> Dict[str, MSExperiment]:
        """Creates an MSExperiment per modification containing all lower-energy spectra
        (listed with their corresponding MS1 spectra to enable searching).

        The result modification identifiers are in amino acid + Unimod accession format,
        windows without a diagnostic ion are grouped into 'unmodified'.

        Windows that contain one of the specified combinations of mods (exact matching)
        are added to the respective combination group. Windows that contain a combination
        that is not listed are added to the groups of the contained single modifications.
        """

        windows_by_mods_combinations = {}
        windows_df_by_mods_combinations = {}
        windows_df_by_mods_single = {}

        for (
            mod_combination,
            spectra_for_mod_df,
        ) in windows_with_detected_ions_df.groupby("letter_and_unimod_format_mod"):
            if mod_combination in mod_combinations_to_search:
                windows_df_by_mods_combinations[mod_combination] = spectra_for_mod_df
            else:
                for mod in mod_combination:
                    if mod not in windows_df_by_mods_single:
                        windows_df_by_mods_single[mod] = [spectra_for_mod_df]
                    else:
                        windows_df_by_mods_single[mod].append(spectra_for_mod_df)

        for mod in windows_df_by_mods_single:
            windows_df_by_mods_single[mod] = pd.concat(windows_df_by_mods_single[mod])

        windows_df_by_mods_combinations = (
            windows_df_by_mods_combinations | windows_df_by_mods_single
        )

        for (
            mod_combination,
            spectra_for_mod_df,
        ) in windows_df_by_mods_combinations.items():
            # TODO: check if one of the sorting columns is actually the index
            spectra_for_mod_df.sort_values(
                ["ms1_spectrum_id", "ms2_spectrum_mz"], inplace=True
            )

            spectra_list_for_mod = np.array([])
            for ms1_spectrum_id, spectra_df in spectra_for_mod_df.groupby(
                ["ms1_spectrum_id"],
                sort=False,  # lexicographical sorting messes up scan window order
            ):
                spectra_list_for_mod = np.append(
                    spectra_list_for_mod,
                    ms1_windows_df.loc[ms1_spectrum_id]["ms1_spectrum"],
                )
                spectra_list_for_mod = np.concatenate(
                    [spectra_list_for_mod, spectra_df["ms2_spectrum"]]
                )

            exp_for_mod = MSExperiment()
            exp_for_mod.setSpectra(list(spectra_list_for_mod))
            windows_by_mods_combinations[mod_combination] = exp_for_mod

        return windows_by_mods_combinations

    def split_windows_by_mods(
        self,
        ms1_windows: MSExperiment,
        ms1_and_lower_energy_windows: MSExperiment,
        ms1_and_higher_energy_windows: MSExperiment,
        detected_ions_df: pd.DataFrame,
        mod_combinations_to_search: List[Tuple[str]],
    ) -> Dict[str, MSExperiment]:
        """Expects the input MSExperiments to be originally extracted from the same
        MSExperiment containing higher- and lower-energy scan windows and the ions
        detected from the higher-energy windows.

        Matches the lower-energy windows to their corresponding higher-energy windows
        to group the lower-energy windows by modifications. Requires the MS1 windows
        as well because the output window list for each modification has to contain the
        MS1 windows before their related MS2 windows in order to be searchable.

        The result modification identifiers are in amino acid + Unimod accession format,
        windows without a diagnostic ion are grouped into 'unmodified'. If no diagnostic
        ions are given, all windows are grouped into 'unmodified'
        """
        self._validate_spectra_mslevel_and_collision_energy(
            ms1_windows, ms1_and_lower_energy_windows, ms1_and_higher_energy_windows
        )

        if len(detected_ions_df) == 0:
            return {"unmodified": ms1_and_lower_energy_windows}

        ms1_windows_df = self._list_ms1_spectra_by_id(ms1_windows)
        num_ms1_windows = len(ms1_windows_df)
        lower_energy_windows_df = self._list_spectra_by_ms1_and_mz(
            ms1_and_lower_energy_windows.getSpectra(), num_ms1_windows
        )
        higher_energy_windows_df = self._list_spectra_by_ms1_and_mz(
            ms1_and_higher_energy_windows.getSpectra(), num_ms1_windows
        )

        self._validate_spectra_identifiers_matching(
            ms1_windows_df, lower_energy_windows_df, higher_energy_windows_df
        )

        higher_energy_windows_detected_ions_df = (
            self._match_higher_energy_spectra_with_detected_ions(
                higher_energy_windows_df, detected_ions_df
            )
        )

        windows_with_detected_ions_df = (
            self._match_lower_and_higher_energy_spectra_with_ions(
                lower_energy_windows_df, higher_energy_windows_detected_ions_df
            )
        )

        windows_by_mods = self._group_windows_by_mods(
            ms1_windows_df, windows_with_detected_ions_df, mod_combinations_to_search
        )

        return windows_by_mods
