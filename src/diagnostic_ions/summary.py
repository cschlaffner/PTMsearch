from pathlib import Path
from typing import FrozenSet, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm


def plot_detected_ions(detected_ions_df: pd.DataFrame, result_path: Path) -> None:
    """Generates a table that lists the immonium ions with their counts.
    The occurences are counted without taking into account whether the occurence is
    in combination with other ions or not."""

    detected_ions_df_by_mod = (
        detected_ions_df[
            [
                "amino_acid",
                "mod_name",
                "letter_and_unimod_format_mod",
                "type",
                "theoretical_mz",
                "spectrum_id",
            ]
        ]
        .rename(columns={"letter_and_unimod_format_mod": "unimod"})
        .groupby(["amino_acid", "mod_name", "unimod", "type", "theoretical_mz"])
        .count()
        .rename(columns={"spectrum_id": "count"})
    )

    fig, ax = plt.subplots()

    # Turn off the figure and only display the table
    # https://stackoverflow.com/a/45936469
    fig.patch.set_visible(False)
    ax.axis("off")
    ax.axis("tight")

    table = ax.table(
        cellText=[
            [str(index_field) for index_field in mod.Index] + [mod.count]
            for mod in detected_ions_df_by_mod.itertuples()
        ],
        colLabels=detected_ions_df_by_mod.index.names
        + [detected_ions_df_by_mod.columns.item()],
        loc="center",
        cellLoc="right",
        colLoc="right",
    )
    # based on https://stackoverflow.com/a/55661458
    cell_height = 1 / len(detected_ions_df_by_mod)
    for _, cell in table.get_celld().items():
        cell.set_height(cell_height)
        cell.set_edgecolor("gray")

    fig.savefig(result_path / "detected_ions_table.svg", bbox_inches="tight")


def get_modification_combinations_with_counts(
    detected_ions_df: pd.DataFrame, return_unimod=False
):
    mod_columns = (
        ["letter_and_unimod_format_mod"]
        if return_unimod
        else ["amino_acid", "mod_name"]
    )
    detected_ions_df_by_window = (
        detected_ions_df[["spectrum_id"] + mod_columns]
        .drop_duplicates()
        .groupby("spectrum_id")
    )

    mods = []
    for _, group in detected_ions_df_by_window:
        mods.append(
            sorted(
                [
                    ",".join(mod)
                    for mod in (group[mod_columns].drop_duplicates().to_numpy())
                ]
            )
        )
    return np.unique(np.array(mods, dtype="object"), return_counts=True)


def plot_detected_ion_combinations(
    detected_ions_df: pd.DataFrame,
    result_path: Path,
    detection_count_percentile: float,
    plot_fractions: List[float] = [0.5, 0.8, 0.9],
    topk: int = 50,
    figsize: Tuple[int, int] = (10, 15),
) -> None:
    combinations, counts = get_modification_combinations_with_counts(
        detected_ions_df, return_unimod=False
    )
    count_sort = np.argsort(counts)
    combinations_sorted = combinations[count_sort]
    combination_counts_sorted = counts[count_sort]
    total_count = combination_counts_sorted.sum()

    if topk is not None:
        combinations_sorted = combinations_sorted[-topk:]
        combination_counts_sorted = combination_counts_sorted[-topk:]

    fig, ax = plt.subplots(figsize=figsize)
    combination_names = [
        (r" $\bf{|}$ ").join(combination) for combination in combinations_sorted
    ]
    num_combinations = len(combination_names)
    ax.bar_label(ax.barh(range(num_combinations), combination_counts_sorted), padding=3)
    ax.set_yticks(range(num_combinations), labels=combination_names)
    ax.set_ylim(0, len(combination_names) + 1)
    ax.set_xlabel("number of windows")
    ax.invert_yaxis()

    plot_fractions.append(detection_count_percentile)
    combination_counts_sorted_reverse = combination_counts_sorted[::-1]
    max_count = combination_counts_sorted_reverse.max()
    color_map = cm.get_cmap("Dark2")

    for j, fraction in enumerate(plot_fractions):
        total_count_fraction = total_count * fraction
        total_count_current = 0

        for i, count in enumerate(combination_counts_sorted_reverse):
            total_count_current += count
            if total_count_current >= total_count_fraction:
                ax.axhline(
                    num_combinations - (i + 1.5), label=fraction, color=color_map(j)
                )
                ax.text(
                    max_count / 2,
                    num_combinations - (i + 1.75),
                    f"Count percentile >= {fraction}",
                    color=color_map(j),
                )
                break

    fig.savefig(
        result_path / f"detected_ions_combinations_top{topk}.svg", bbox_inches="tight"
    )


def get_detected_modifications_with_combinations(
    detected_ions_df: pd.DataFrame,
    detection_count_percentile: float,
    detection_count_min: int,
    num_additional_modifications: int,
) -> Tuple[List[str], FrozenSet[str]]:

    combinations, counts = get_modification_combinations_with_counts(
        detected_ions_df, return_unimod=True
    )
    count_sort = np.argsort(counts)
    combinations_sorted = combinations[count_sort][::-1]
    combination_counts_sorted = counts[count_sort][::-1]
    total_count = combination_counts_sorted.sum()

    combinations_sorted_min = combinations_sorted[
        combination_counts_sorted >= detection_count_min
    ]
    combination_counts_sorted_min = combination_counts_sorted[
        combination_counts_sorted >= detection_count_min
    ]

    total_count_fraction = total_count * detection_count_percentile
    total_count_current = 0

    percentile_index = len(combination_counts_sorted_min)
    for i, count in enumerate(combination_counts_sorted_min):
        total_count_current += count
        if total_count_current >= total_count_fraction:
            percentile_index = i
            break

    combinations_percentile = combinations_sorted_min[: percentile_index + 1]
    # DIA-NN can only handle max. 5 modifications at once
    # and single mods should not be handled as combination
    combinations_sets = [
        frozenset(combination)
        for combination in combinations_percentile
        if len(combination) + num_additional_modifications <= 5 and len(combination) > 1
    ]

    modifications = np.unique(np.concatenate(combinations_sorted))
    modifications_min = []

    discarded_combinations_with_counts = [
        (combination, count)
        for combination, count in zip(combinations_sorted, combination_counts_sorted)
        if frozenset(combination) not in combinations_sets
    ]
    for mod in modifications:
        count_discarded_combinations_mod = np.sum(
            [
                count
                for combination, count in discarded_combinations_with_counts
                if mod in combination
            ]
        )

        if count_discarded_combinations_mod >= np.maximum(detection_count_min, 1):
            modifications_min.append(mod)

    return modifications_min, combinations_sets
