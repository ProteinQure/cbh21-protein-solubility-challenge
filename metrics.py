import numpy as np
import pandas as pd
import sys
from scipy.stats import spearmanr, pearsonr

from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_percentage_error
from sklearn.metrics import mean_squared_log_error
"""
def compute_metrics(predictions_file: str, target_files: dict[str, str]) -> dict[str, Any]:
    
    #Computes the test performance metrics on the three test data sets in
    #`target_files`.
    

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
"""
y_true = pd.read_csv('True.csv', index_col=False)
y_pred = pd.read_csv('predictions.csv', index_col=False)

y_true = y_true.iloc[:,-1].tolist()
y_pred = y_pred.iloc[:,-1].tolist()

y_true = [float(i) for i in y_true]

print(y_true)
print(y_pred)

print(mean_absolute_percentage_error(y_true, y_pred))
print(explained_variance_score(y_true, y_pred))
print(mean_absolute_error(y_true, y_pred))

rmse = mean_squared_error(y_true, y_pred, squared=False)
spr = spearmanr(a=y_true,b=y_pred)
pec = pearsonr(x=y_true,y=y_pred)

print(rmse)
print(spr)
print(pec)

#1123.7541
#32 = old RMSE