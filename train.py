from feature_construction import compute_features
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import glob
import csv
import pickle
import freesasa


# train the random forrest with training data

filenames = glob.glob("data/training/crystal_structs/*.pdb")

features = compute_features(filenames)

# features = np.load('features.npy')
# features = pd.read_csv('features.csv')
features = features.sort_values('protIDs').reset_index(drop=True)

solubility = pd.read_csv("data/training/crystal_structs/solubility_values.csv")
solubility = solubility.sort_values('protein').reset_index(drop=True)

# merge solubility and features (not in the same order?)
# merged = pd.merge(features, solubility, left_on='protIDs', right_on='protein')

# save the proteins 
proteins = solubility['protein'].values

X = features.iloc[:, 1:].to_numpy()
y = solubility.iloc[:, 1:].to_numpy()

# TrainX, TestX, TrainY, TestY, proteinX, proteinY = solubility

XTrain, XTest, YTrain, YTest = train_test_split(
    X, y, test_size=0.33, random_state=42)

rfr = RandomForestRegressor(n_estimators=500)
rfr.fit(XTrain, YTrain)
rfr_prediction = rfr.predict(XTest)

from sklearn import metrics

print(rfr.feature_importances_)

print('Mean Absolute Error (MAE):', metrics.mean_absolute_error(YTest, rfr_prediction))
print('Mean Squared Error (MSE):', metrics.mean_squared_error(YTest, rfr_prediction))
print('Root Mean Squared Error (RMSE):', np.sqrt(metrics.mean_squared_error(YTest, rfr_prediction)))
mape = np.mean(np.abs((YTest - rfr_prediction) / np.abs(YTest)))
print('Mean Absolute Percentage Error (MAPE):', round(mape * 100, 2))
print('Accuracy:', round(100*(1 - mape), 2))



# save model

pkl_filename = "pickle_model.pkl"
with open(pkl_filename, 'wb') as file:
    pickle.dump(rfr, file)
