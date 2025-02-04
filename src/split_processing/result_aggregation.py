import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MzMLFile


class ResultAggregation:
    """A class to aggregate results from splits and calculate aggregated q-values."""

    def __init__(self, fdr_threshold: float, normalize_cscores: bool = False):
        self.fdr_threshold = fdr_threshold
        self.normalize_cscores = normalize_cscores
        self.c_score_column = "CScore_normalized" if normalize_cscores else "CScore"

    def _fix_decoy_report(self, df):
        """The decoy precursor rows in the DIA-NN report are misaligned and therefore
        columns do not match, so the relevant column values must be shifted to
        the right columns after loading. Further, decoys do not have their window
        of detection stored, so the respective final report columns have to get
        placeholders.
        """
        df.loc[:, "Q.Value"] = pd.to_numeric(df["Precursor.Id"])
        df.loc[:, "Precursor.Id"] = df["Modified.Sequence"]
        df.loc[:, "CScore"] = df["RT.Start"]

        df.insert(len(df.columns), "remapped_id", np.repeat(np.nan, len(df)))
        df.insert(len(df.columns), "original_id", np.repeat(np.nan, len(df)))

    def _get_id_number(self, id_string):
        return int(re.findall("[0-9]+", id_string)[-1])

    def _get_df_with_original_id_mapped(self, df, exp, mapping_df):
        """Recovers the ID of the renamed spectrum from the search mzML and
        the original ID that the spectrum had in the original dataset without
        splitting from the DIA-NN report column that lists the index of the MS2 scan
        (considering only MS2 scans for the internal indexing apparently, see
        https://github.com/vdemichev/DiaNN/discussions/985)."""
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

    def _qvalues_for_cscores(
        self, cscores, targets_cscores, decoys_cscores, batch_size=1000
    ):
        """Calculates q-values as local FDR for the given c-scores.
        Formula follows the official DIA-NN paper (Demichev et al. (2017)).
        Calculation is done in batches to make the procedure a bit faster."""
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

    def _qvalues_for_targets_and_decoys(
        self, targets_df: pd.DataFrame, decoys_df: pd.DataFrame
    ):
        """Wraps q-value calculation for targets and decoys."""
        targets_cscores = targets_df[self.c_score_column].to_numpy()
        decoys_cscores = decoys_df[self.c_score_column].to_numpy()
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
        """Min-max c-score normalization to 0-1 range. DIA-NN sometimes
        produces c-scores larger than 1 for a split, which
        can skew the overall c-score distributions. This can be mitigated
        by normalization."""
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

        split_targets_df.insert(
            len(split_targets_df.columns),
            self.c_score_column,
            target_cscores_normalized,
        )
        split_decoys_df.insert(
            len(split_decoys_df.columns), self.c_score_column, decoy_cscores_normalized
        )

    def _de_duplicate_precursors(self, targets_df: pd.DataFrame):
        """Removes duplicate precursors, i.e., precursors with same ID found
        in the same window, while keeping the one with the maximum c-score."""
        # based on https://stackoverflow.com/a/45527762
        return (
            targets_df.sort_values(self.c_score_column, ascending=False)
            .groupby(by=["Precursor.Id", "original_id"], as_index=False)
            .head(1)
            .reset_index(drop=True)
        )

    def _get_mods_unmods_all_from_splits(self, file_paths_by_mods):
        """Obtains precursors (targets and decoys) by splits:
        all modified splits together, the unmodified split
        separately, and all splits together."""

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
            self._get_df_with_original_id_mapped(targets, exp, mapping_df)

            if self.normalize_cscores:
                self._normalize_cscores(targets, decoys)

            splits_results[mods] = {"targets": targets, "decoys": decoys}

        mods_splits = {
            mods: results
            for mods, results in splits_results.items()
            if mods != "unmodified"
        }

        if len(mods_splits) > 0:
            mods_targets = pd.concat(
                [mod_result["targets"] for mod_result in mods_splits.values()],
                ignore_index=True,
            )

            mods_decoys = pd.concat(
                [mod_result["decoys"] for mod_result in mods_splits.values()],
                ignore_index=True,
            )
        else:
            mods_targets = pd.DataFrame()
            mods_decoys = pd.DataFrame()

        if "unmodified" in splits_results:
            unmods_targets = splits_results["unmodified"]["targets"]
            unmods_decoys = splits_results["unmodified"]["decoys"]
        else:
            unmods_targets = pd.DataFrame()
            unmods_decoys = pd.DataFrame()

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
        """Density plots for the c-score distributions. Density plotting
        only succeeds when the data has at least 2 elements."""
        if len(targets_df) <= 1 or len(decoys_df) <= 1:
            return

        if normalized:
            decoy_cscores = decoys_df[self.c_score_column]
            target_cscores = targets_df[self.c_score_column]
            cscore_label = "cscore_normalized"
        else:
            decoy_cscores = decoys_df["CScore"]
            target_cscores = targets_df["CScore"]
            cscore_label = "cscore_unnormalized"

        plt.figure()
        # to have unified plots for c-scores in the expected normal range
        if (
            decoy_cscores.max() < 1.5
            and decoy_cscores.min() >= 0
            and target_cscores.max() < 1.5
            and target_cscores.min() >= 0
        ):
            plt.xlim((-0.25, 1.5))
            plt.xticks(np.arange(0, 1.26, 0.25))

        decoy_cscores.plot.kde(label="Decoys")
        target_cscores.plot.kde(label="Targets")
        plt.xlabel("CScore")
        plt.legend()

        plt.savefig(
            result_path / f"{cscore_label}_densities_{fig_name}.png",
            bbox_inches="tight",
        )

    def aggregate_results(self, result_path, file_paths_by_mods):
        """Aggregates the DIA-NN reports from the different splits
        and returns a full report (no FDR filtering, containing targets
        and decoys) and a filtered report (only targets that pass FDR
        filtering and duplicates removed). Includes plotting of the
        c-score distributions."""
        (
            mods_targets,
            mods_decoys,
            unmods_targets,
            unmods_decoys,
            all_targets,
            all_decoys,
        ) = self._get_mods_unmods_all_from_splits(file_paths_by_mods)

        assert len(all_targets) > 1 and len(all_decoys) > 1, (
            "At most 1 targets and/or at most 1 decoys were listed in the "
            "results even before FDR filtering. Looks like something went "
            "wrong regarding configurations or data for this run."
        )

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

        if self.normalize_cscores:
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

        self._qvalues_for_targets_and_decoys(all_targets, all_decoys)

        report_aggregated = pd.concat(
            [all_targets, all_decoys],
            ignore_index=True,
        )

        cscore_threshold = all_targets[
            all_targets["q_value_aggregated"] <= self.fdr_threshold
        ][self.c_score_column].min()

        all_targets_no_duplicates = self._de_duplicate_precursors(all_targets)

        report_fdr_filtered = all_targets_no_duplicates[
            all_targets_no_duplicates[self.c_score_column] >= cscore_threshold
        ]

        return report_aggregated, report_fdr_filtered
