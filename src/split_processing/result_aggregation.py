# TODO This will contain the result merging, regarding FDR and stuff, and (probably) protein inference

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MzMLFile


class ResultAggregation:
    def __init__(self, fdr_threshold: float):
        self.fdr_threshold = fdr_threshold

    def _fix_decoy_report(self, df):
        df.loc[:, "Q.Value"] = pd.to_numeric(df["Precursor.Id"])
        df.loc[:, "Precursor.Id"] = df["Modified.Sequence"]
        df.loc[:, "CScore"] = df["RT.Start"]

        df.insert(len(df.columns), "remapped_id", np.repeat(np.nan, len(df)))
        df.insert(len(df.columns), "original_id", np.repeat(np.nan, len(df)))
        df.insert(
            len(df.columns),
            "higher_energy_id",
            np.repeat(np.nan, len(df)),
        )

    def _get_id_number(self, id_string):
        return int(re.findall("[0-9]+", id_string)[-1])

    def _get_higher_energy_scan_id(self, id_string):
        id_number = self._get_id_number(id_string)
        higher_energy_id_number = id_number + 36
        return f"controllerType=0 controllerNumber=1 scan={higher_energy_id_number}"

    def _get_df_with_original_and_higher_energy_id_mapped(self, df, exp, mapping_df):
        ms2_spectra = np.array(
            [spectrum for spectrum in exp.getSpectra() if spectrum.getMSLevel() == 2]
        )
        ms2_spectra_in_df = ms2_spectra[df["MS2.Scan"].to_numpy().astype(int)]
        df.insert(
            len(df.columns),
            "remapped_id",
            [s.getNativeID() for s in ms2_spectra_in_df],
        )
        df.insert(
            len(df.columns),
            "original_id",
            [
                mapping_df.loc[remapped_id]["original_id"]
                for remapped_id in df["remapped_id"]
            ],
        )
        higher_energy_id_number = np.array(
            [
                self._get_id_number(self._get_higher_energy_scan_id(id))
                for id in df["original_id"]
            ]
        )
        df.insert(len(df.columns), "higher_energy_id", higher_energy_id_number)

    def _qvalues_for_cscores(
        self, cscores, targets_cscores, decoys_cscores, batch_size=1000
    ):
        cscores = cscores[:, None]
        targets_cscores = targets_cscores[None, :]
        decoys_cscores = decoys_cscores[None, :]

        split_borders = np.arange(0, len(cscores) + batch_size, batch_size)
        q_values = []

        for i in range(len(split_borders) - 1):
            current_cscores = cscores[split_borders[i] : split_borders[i + 1]]
            targets_passing = targets_cscores > current_cscores
            decoys_passing = decoys_cscores > current_cscores
            q_values.append(decoys_passing.sum(axis=1) / targets_passing.sum(axis=1))

        return np.concatenate(q_values)

    def _compute_qvalues(self, targets_df: pd.DataFrame, decoys_df: pd.DataFrame):
        targets_cscores = targets_df["CScore"].to_numpy()
        decoys_cscores = decoys_df["CScore"].to_numpy()
        qvalues_targets = self._qvalues_for_cscores(
            targets_cscores, targets_cscores, decoys_cscores
        )
        qvalues_decoys = self._qvalues_for_cscores(
            decoys_cscores, targets_cscores, decoys_cscores
        )

        targets_df.insert(
            len(targets_df.columns), "q_value_aggregated", qvalues_targets
        )
        decoys_df.insert(len(decoys_df.columns), "q_value_aggregated", qvalues_decoys)

    def _normalize_cscores(
        self, split_targets_df: pd.DataFrame, split_decoys_df: pd.DataFrame
    ):
        target_cscores = split_targets_df["CScore"].to_numpy()
        decoy_cscores = split_decoys_df["CScore"].to_numpy()
        split_cscores = np.concatenate([target_cscores, decoy_cscores])
        scores_min = split_cscores.min()
        scores_max = split_cscores.max()

        # Min-max normalization to 0-1 range
        target_cscores_normalized = (target_cscores - scores_min) / (
            scores_max - scores_min
        )
        decoy_cscores_normalized = (decoy_cscores - scores_min) / (
            scores_max - scores_min
        )

        split_targets_df.rename(columns={"CScore": "CScore_unnormalized"}, inplace=True)
        split_targets_df.insert(
            len(split_targets_df.columns), "CScore", target_cscores_normalized
        )
        split_decoys_df.rename(columns={"CScore": "CScore_unnormalized"}, inplace=True)
        split_decoys_df.insert(
            len(split_decoys_df.columns), "CScore", decoy_cscores_normalized
        )

    def _aggregate_duplicate_precursors(self, targets_df: pd.DataFrame):
        # based on https://stackoverflow.com/a/45527762
        return (
            targets_df.sort_values("CScore", ascending=False)
            .groupby(by=["Precursor.Id", "original_id"], as_index=False)
            .head(1)
            .reset_index(drop=True)
        )

    def get_mods_unmods_all_from_splits(self, file_paths_by_mods, normalize_cscores):
        splits_results = {}
        for mods in file_paths_by_mods:
            file_paths = file_paths_by_mods[mods]

            dia_nn_report_df = pd.read_csv(
                file_paths["dia_nn_report_path"], delimiter="\t"
            )
            mapping_df = pd.read_csv(
                file_paths["spectrum_id_mapping_path"], index_col="renamed_id"
            )

            exp = MSExperiment()
            MzMLFile().load(str(file_paths["mzml_path"]), exp)

            # to account for mismatching columns in result tsv for decoys
            decoys = dia_nn_report_df[dia_nn_report_df["Q.Value"].isna()]
            self._fix_decoy_report(decoys)

            targets = dia_nn_report_df[~dia_nn_report_df["Q.Value"].isna()]
            self._get_df_with_original_and_higher_energy_id_mapped(
                targets, exp, mapping_df
            )

            if normalize_cscores:
                self._normalize_cscores(targets, decoys)

            splits_results[mods] = {"targets": targets, "decoys": decoys}

        mods_splits = {
            mods: results
            for mods, results in splits_results.items()
            if mods != "unmodified"
        }

        mods_targets = pd.concat(
            [mod_result["targets"] for mod_result in mods_splits.values()],
            ignore_index=True,
        )
        mods_decoys = pd.concat(
            [mod_result["decoys"] for mod_result in mods_splits.values()],
            ignore_index=True,
        )

        unmods_targets = splits_results["unmodified"]["targets"]
        unmods_decoys = splits_results["unmodified"]["decoys"]

        all_targets = pd.concat([mods_targets, unmods_targets], ignore_index=True)
        all_decoys = pd.concat([mods_decoys, unmods_decoys], ignore_index=True)

        return (
            mods_targets,
            mods_decoys,
            unmods_targets,
            unmods_decoys,
            all_targets,
            all_decoys,
        )

    def plot_densities(
        self, result_path, targets_df, decoys_df, fig_name, normalized=False
    ):
        if normalized:
            decoy_cscores = decoys_df["CScore"]
            target_cscores = targets_df["CScore"]
            cscore_label = "cscore_normalized"
        else:
            decoy_cscores = decoys_df["CScore_unnormalized"]
            target_cscores = targets_df["CScore_unnormalized"]
            cscore_label = "cscore_unnormalized"

        plt.figure()
        decoy_cscores.plot.kde(label="decoys")
        target_cscores.plot.kde(label="targets")
        plt.xlabel("CScore")
        plt.legend()

        plt.savefig(
            result_path / f"{cscore_label}_densities_{fig_name}.png",
            bbox_inches="tight",
        )

    def aggregate_results(
        self, result_path, file_paths_by_mods, normalize_cscores=False
    ):
        (
            mods_targets,
            mods_decoys,
            unmods_targets,
            unmods_decoys,
            all_targets,
            all_decoys,
        ) = self.get_mods_unmods_all_from_splits(file_paths_by_mods, normalize_cscores)

        self.plot_densities(
            result_path,
            mods_targets,
            mods_decoys,
            "modification_splits",
            normalized=False,
        )
        self.plot_densities(
            result_path,
            unmods_targets,
            unmods_decoys,
            "unmodified_split",
            normalized=False,
        )
        self.plot_densities(
            result_path, all_targets, all_decoys, "all_splits", normalized=False
        )

        if normalize_cscores:
            self.plot_densities(
                result_path,
                mods_targets,
                mods_decoys,
                "modification_splits",
                normalized=True,
            )
            self.plot_densities(
                result_path,
                unmods_targets,
                unmods_decoys,
                "unmodified_split",
                normalized=True,
            )
            self.plot_densities(
                result_path, all_targets, all_decoys, "all_splits", normalized=True
            )

        self._compute_qvalues(all_targets, all_decoys)

        report_aggregated = pd.concat(
            [all_targets, all_decoys],
            ignore_index=True,
        )

        all_targets_no_duplicates = self._aggregate_duplicate_precursors(all_targets)

        report_fdr_filtered = all_targets_no_duplicates[
            all_targets_no_duplicates["q_value_aggregated"] <= self.fdr_threshold
        ]

        return report_aggregated, report_fdr_filtered
