from typing import Dict, List

import numpy as np
import pandas as pd
from pyopenms import MSSpectrum

from src.mzml_processing.utils import check_ms1_spectrum


class ScanWindowSplitting:
    def list_spectra_by_ms1_and_rt(
        self,
        ms1_and_ms2_spectra: List[MSSpectrum],
    ) -> pd.DataFrame:
        """Creates a DataFrame where every MS2 spectrum is identified by its retention time and corresponding
        MS1 spectrum. This is required to match higher- and lower-energy spectra."""
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

    def split_windows_by_mods(
        self,
        ms1_and_lower_energy_windows: List[MSSpectrum],
        ms1_and_higher_energy_windows: List[MSSpectrum],
        detected_ions_df: pd.DataFrame,
    ) -> Dict[str, List[MSSpectrum]]:

        lower_energy_windows_df = self.list_spectra_by_ms1_and_rt(
            ms1_and_lower_energy_windows
        )

        lower_energy_windows_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_rt"], inplace=True
        )

        higher_energy_windows_df = self.list_spectra_by_ms1_and_rt(
            ms1_and_higher_energy_windows
        )

        higher_energy_windows_df.set_index("ms2_spectrum_id", inplace=True)

        # Here, we do not need any additional information (e.g. the ion type) anymore
        detected_ions_df = detected_ions_df[
            ["spectrum_id", "amino_acid", "mod_name"]
        ].drop_duplicates()

        # TODO: convert NaN entries to "unmodified"
        higher_energy_windows_detected_ions_df = detected_ions_df.join(
            higher_energy_windows_df, on="spectrum_id", how="right"
        )

        # TODO: reset index beforehand?
        higher_energy_windows_detected_ions_df.set_index(
            ["ms1_spectrum_id", "ms2_spectrum_rt"], inplace=True
        )
        windows_with_detected_mods_df = higher_energy_windows_detected_ions_df.join(
            lower_energy_windows_df,
            on=["ms1_spectrum_id", "ms2_spectrum_rt"],
            how="inner",  # although that shouldn't matter here
        )

        window_list_by_mods = {}

        for mod_name, spectra_for_mod_df in windows_with_detected_mods_df.groupby(
            ["amino_acid", "mod_name"]
        ):
            spectra_for_mod_list = []
            for ms1_spectrum, spectra_df in spectra_for_mod_df.groupby(
                ["ms1_spectrum_id", "ms1_spectrum"]
            ):
                spectra_for_mod_list.append(ms1_spectrum[1])
                spectra_for_mod_list = np.concatenate(
                    spectra_for_mod_list, spectra_df["ms2_spectrum"]
                )

            # TODO: check if it works with tuple key
            window_list_by_mods[mod_name] = spectra_for_mod_list
        return window_list_by_mods
