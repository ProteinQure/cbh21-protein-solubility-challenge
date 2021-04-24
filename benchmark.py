"""
Script used for the automatic benchmarking of the challenge.
"""

from __future__ import annotations
from typing import Any

import pandas as pd
import argparse
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr


def compute_metrics(predictions_file: str, target_files: list) -> dict[str, Any]:
    """
    Computes the test performance metrics on the test data sets in
    `target_files`.
    """

    predictions_df = pd.read_csv(predictions_file)
    predictions_df.columns = ["protein", "prediction"]
    targets_dfs = [pd.read_csv(target_file) for target_file in target_files]
    targets_df = pd.concat(targets_dfs)
    targets_df.columns = ["protein", "true_value"]

    # merge the two dataframes
    df = pd.merge(predictions_df, targets_df, on="protein")

    true_values = df["true_value"].to_numpy()
    predictions = df["prediction"].to_numpy()
    mse = mean_squared_error(y_true=true_values, y_pred=predictions)
    correlation, pvalue = spearmanr(a=true_values, b=predictions)

    return {
        "mean_squared_error": mse,
        "spearman_rank": {
            "correlation": correlation,
            "pvalue": pvalue,
        },
        "num_points_used_for_metrics": f"{df.shape[0]}/{targets_df.shape[0]}",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", default="predictions.csv")
    parser.add_argument("--yeast-crystal", default="data/test/yeast_crystal_structs/solubility_values.csv")
    parser.add_argument("--yeast-modelled", default="data/test/yeast_modelled_structs/solubility_values.csv")
    parser.add_argument("--ecoli-modelled", default="data/test/ecoli_modelled_structs/solubility_values.csv")
    args = parser.parse_args()
    print(
        "Test set performance: ",
        compute_metrics(args.predictions, [args.yeast_crystal, args.yeast_modelled, args.ecoli_modelled]),
    )
