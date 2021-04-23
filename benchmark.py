"""
Script used for the automatic benchmarking of the challenge.
"""

from __future__ import annotations
import argparse
import pprint
from typing import Any

import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr


def compute_metrics(predictions_file: str, target_files: dict[str, str]) -> dict[str, Any]:
    """
    Computes the test performance metrics on the three test data sets in
    `target_files`.
    """

    predictions_df = pd.read_csv(predictions_file)
    predictions_df.columns = ["protein", "prediction"]

    metrics = {}
    for test_set_name, target_file in target_files.items():

        # this allows users to run it with the `yeast_crystal_structs` test set only
        if target_file is None:
            continue

        targets_df = pd.read_csv(target_file)
        targets_df.columns = ["protein", "true_value"]

        # merge the two dataframes (gets rid of the rows we don't need)
        df = pd.merge(predictions_df, targets_df, on="protein")

        true_values = df["true_value"].to_numpy()
        predictions = df["prediction"].to_numpy()
        rmse = mean_squared_error(y_true=true_values, y_pred=predictions, squared=False)
        spearman_correlation, spearman_pvalue = spearmanr(a=true_values, b=predictions)
        pearson_correlation, pearson_pvalue = pearsonr(x=true_values, y=predictions)

        metrics.update(
            {
                test_set_name: {
                    "rmse": rmse,
                    "spearman_rank": {
                        "correlation": spearman_correlation,
                        "pvalue": spearman_pvalue,
                    },
                    "pearson": {
                        "correlation": pearson_correlation,
                        "pvalue": pearson_pvalue,
                    },
                }
            }
        )

    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", default="predictions.csv")
    parser.add_argument("--yeast-crystal", default="data/test/yeast_crystal_structs/solubility_values.csv")
    parser.add_argument("--yeast-modelled", default=None)
    parser.add_argument("--ecoli-modelled", default=None)
    args = parser.parse_args()

    pp = pprint.PrettyPrinter(indent=4, width=120)
    print("TEST SET PERFORMANCES")
    print("=====================")
    pp.pprint(
        compute_metrics(
            args.predictions,
            {
                "yeast_crystal_structs": args.yeast_crystal,
                "yeast_modelled_structs": args.yeast_modelled,
                "ecoli_modelled_structs": args.ecoli_modelled,
            },
        ),
    )
