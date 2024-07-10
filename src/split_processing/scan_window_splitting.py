from typing import Dict, List

import numpy as np
import pandas as pd
import pyopenms as oms
from pyopenms import MSExperiment, MSSpectrum

from src.mzml_processing.utils import check_ms1_spectrum


class ScanWindowSplitting:

    def __init__(self) -> None:
        self.modifications_db = oms.ModificationsDB()
        self.residue_db = oms.ResidueDB()

    def list_spectra_by_ms1_and_rt(
        self,
        ms1_and_ms2_spectra: List[MSSpectrum],
    ) -> pd.DataFrame:
        """Creates a DataFrame where every MS2 spectrum is identified by its retention time and corresponding
        MS1 spectrum ID. This is required to match higher- and lower-energy spectra."""
        spectra_empty_df = pd.DataFrame(
            columns=[
                "ms1_spectrum_id",
                "ms2_spectrum_rt",
                "ms2_spectrum_id",
                "ms2_spectrum",
            ]
        )
        spectra_dfs = [spectra_empty_df]
        assert check_ms1_spectrum(
            ms1_and_ms2_spectra[0]
        ), "Scan window list should start with an MS1 spectrum."
        current_ms1_spectrum = None

        for spectrum in ms1_and_ms2_spectra:
            if check_ms1_spectrum(spectrum):
                current_ms1_spectrum = spectrum
            else:
                spectra_dfs.append(
                    pd.DataFrame(
                        {
                            "ms1_spectrum_id": current_ms1_spectrum.getNativeID(),
                            "ms2_spectrum_rt": spectrum.getRT(),
                            "ms2_spectrum_id": spectrum.getNativeID(),
                            "ms2_spectrum": spectrum,
                        }
                    )
                )

        spectra_df = pd.concat(spectra_dfs, ignore_index=True)
        return spectra_df

    def list_ms1_spectra_by_id(self, ms1_windows: MSExperiment) -> pd.DataFrame:
        """Helper function because normal MSExperiment does not allow retrieving spectra by native ID."""
        spectra = ms1_windows.getSpectra()
        spectra_df = pd.DataFrame(
            columns=["ms1_spectrum_id", "ms1_spectrum"], index=range(len(spectra))
        )
        for i, spectrum in enumerate(spectra):
            spectra_df.loc[i] = {
                "ms1_spectrum_id": spectrum.getNativeID(),
                "ms1_spectrum": spectrum,
            }
        spectra_df.set_index("ms1_spectrum_id", inplace=True)
        return spectra_df

    def _validate_spectra_identifiers_matching(
        self,
        ms1_windows_df: pd.DataFrame,
        lower_energy_windows_df: pd.DataFrame,
        higher_energy_windows_df: pd.DataFrame,
    ) -> None:
        assert lower_energy_windows_df[["ms1_spectrum_id", "ms2_spectrum_rt"]].equals(
            higher_energy_windows_df[["ms1_spectrum_id", "ms2_spectrum_rt"]]
        )

        assert ms1_windows_df["ms1_spectrum_id"].equals(
            lower_energy_windows_df["ms1_spectrum_id"].drop_duplicates()
        )

    def match_higher_energy_spectra_with_detected_ions(
        self, higher_energy_windows_df: pd.DataFrame, detected_ions_df: pd.DataFrame
    ) -> pd.DataFrame:
        """The data about the detected ions contains only the ID of the higher-energy spectrum.
        To match the ions with the lower-energy spectra, more higher-energy window information
        (retention time and MS1 spectrum ID) is needed.
        """

        higher_energy_windows_df.set_index("ms2_spectrum_id", inplace=True)

        # Here, we do not need any additional information (e.g. the ion type) anymore
        detected_ions_df = detected_ions_df[
            ["spectrum_id", "amino_acid", "mod_name"]
        ].drop_duplicates()

        higher_energy_windows_detected_ions_df = detected_ions_df.join(
            higher_energy_windows_df, on="spectrum_id", how="right"
        )
        higher_energy_windows_detected_ions_df.loc[
            higher_energy_windows_detected_ions_df["mod_name"].isna(), "mod_name"
        ] = "unmodified"

        # TODO: reset index beforehand?
        return higher_energy_windows_detected_ions_df

    def match_lower_and_higher_energy_spectra_with_ions(
        self,
        lower_energy_windows_df: pd.DataFrame,
        higher_energy_windows_detected_ions_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Matches the lower-energy spectra with the information about the detected ions
        containing also retention time and MS1 spectrum ID for each higher-energy window
        in order to perform the matching.
        """

        # re-index to improve performance during join
        lower_energy_windows_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_rt"], inplace=True
        )
        # TODO: reset index beforehand?
        higher_energy_windows_detected_ions_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_rt"], inplace=True
        )

        windows_with_detected_ions_df = higher_energy_windows_detected_ions_df.join(
            lower_energy_windows_df,
            on=["ms1_spectrum_id", "ms2_spectrum_rt"],
            how="inner",  # although that shouldn't matter here
        )
        return windows_with_detected_ions_df

    def group_windows_by_mods(
        self, ms1_windows_df: pd.DataFrame, windows_with_detected_ions_df: pd.DataFrame
    ) -> Dict[str, MSExperiment]:
        """Creates an MSExperiment per modification containing all lower-energy spectra
        (listed with their corresponding MS1 spectra to enable searching).

        The result modification identifiers are in amino acid + Unimod accession format, windows without a diagnostic ion
        are grouped into 'unmodified'.
        """

        windows_by_mods = {}

        for mod, spectra_for_mod_df in windows_with_detected_ions_df.groupby(
            ["amino_acid", "mod_name"]
        ):
            spectra_list_for_mod = np.array([])
            for ms1_spectrum_id, spectra_df in spectra_for_mod_df.groupby(
                ["ms1_spectrum_id"]
            ):
                spectra_list_for_mod = np.append(
                    spectra_list_for_mod,
                    ms1_windows_df.loc[ms1_spectrum_id]["ms1_spectrum"],
                )
                spectra_list_for_mod = np.concatenate(
                    [spectra_list_for_mod, spectra_df["ms2_spectrum"]]
                )

            exp_for_mod = MSExperiment()
            exp_for_mod.setSpectra(spectra_list_for_mod)

            amino_acid, mod_name = mod
            if mod_name == "unmodified":
                windows_by_mods[mod_name] = exp_for_mod
            else:
                unimod_accession = self.modifications_db.getModification(
                    mod_name
                ).getUniModAccession()
                amino_acid_letter = self.residue_db.getResidue(
                    amino_acid
                ).getOneLetterCode()
                windows_by_mods[amino_acid_letter + unimod_accession] = exp_for_mod

        return windows_by_mods

    def split_windows_by_mods(
        self,
        ms1_windows: MSExperiment,
        ms1_and_lower_energy_windows: MSExperiment,
        ms1_and_higher_energy_windows: MSExperiment,
        detected_ions_df: pd.DataFrame,
    ) -> Dict[str, MSExperiment]:
        """Expects the input MSExperiments to be originally extracted from the same MSExperiment containing higher-
        and lower-energy scan windows and the ions detected from the higher-energy windows.

        Matches the lower-energy windows to their corresponding higher-energy windows to group the lower-energy windows
        by modifications. Requires the MS1 windows as well because the output window list for each modification has to
        contain the MS1 windows before their related MS2 windows in order to be searchable.

        The result modification identifiers are in amino acid + Unimod accession format, windows without a diagnostic ion
        are grouped into 'unmodified'.
        """

        ms1_windows_df = self.list_ms1_spectra_by_id(ms1_windows)
        lower_energy_windows_df = self.list_spectra_by_ms1_and_rt(
            ms1_and_lower_energy_windows.getSpectra()
        )
        higher_energy_windows_df = self.list_spectra_by_ms1_and_rt(
            ms1_and_higher_energy_windows.getSpectra()
        )

        self._validate_spectra_identifiers_matching(
            ms1_windows_df, lower_energy_windows_df, higher_energy_windows_df
        )

        higher_energy_windows_detected_ions_df = (
            self.match_higher_energy_spectra_with_detected_ions(
                higher_energy_windows_df, detected_ions_df
            )
        )

        windows_with_detected_ions_df = (
            self.match_lower_and_higher_energy_spectra_with_ions(
                lower_energy_windows_df, higher_energy_windows_detected_ions_df
            )
        )

        windows_by_mods = self.group_windows_by_mods(
            ms1_windows_df, windows_with_detected_ions_df
        )

        return windows_by_mods
