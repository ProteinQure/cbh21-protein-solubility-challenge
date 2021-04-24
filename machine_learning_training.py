import numpy as np
import pandas as pd
import sys

from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn import datasets, ensemble
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.model_selection import KFold
from sklearn import metrics
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn import metrics
import time
import warnings
import pickle


np.set_printoptions(threshold=np.inf)

data_train_features = pd.read_csv("training_data.csv")
#data_train_features['0'] = data_train_features.index
#del data_train_features['0']

data_train_outcome = pd.read_csv("solubility_train.csv")

data_test_features = pd.read_csv("yeast_crystal.csv")
#data_test_features['0'] = data_train_features.index
#del data_test_features['0']

data_test_outcome = pd.read_csv("solubility_yeast_crystal.csv")

X_train = data_train_features
y_train = np.asarray(data_train_outcome.iloc[:, -1])

X_test = data_test_features
y_test = np.asarray(data_test_outcome.iloc[:,-1])


#'GroupA', 'GroupB',, 'PrimOut_di'
'''discrete_transformer = Pipeline(steps=[
    ('imputer', SimpleImputer(missing_values=12345, strategy='most_frequent')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))])'''

continuous_data = data_train_features.columns.tolist()[:-1]

#encoder = preprocessing.StandardScaler()
encoder = preprocessing.KBinsDiscretizer(n_bins=10, encode='onehot', strategy='kmeans')
continuous_transformer = Pipeline(steps=[
    #('imputer', SimpleImputer(missing_values=12345, strategy='mean')),
    ('discretizer', encoder)])

#decides which transformer to use based on the whether the column is in the discrete or contiuous list
preprocessor = ColumnTransformer(
    transformers=[
       ('continuous', continuous_transformer, continuous_data)])

le = preprocessing.LabelEncoder()

class_weight = {}
regressor = RandomForestRegressor(n_estimators=200,
                                  oob_score=True,
                                  bootstrap=True,
                                  random_state=33)#, class_weight=class_weight)

classifier = RandomForestClassifier(n_estimators=200, oob_score=True, bootstrap=True, class_weight=class_weight)

selector = RFE(classifier, n_features_to_select = 10, step=0.10)

#pipeline for the GridSearchCV to use
pipeline = Pipeline(steps=[("preprocessor", preprocessor),
                            ("RFE", selector),
                            ("regressor", regressor)])

#parameter grid for the GridSearchCV to use
param_grid = dict(RFE__n_features_to_select=[3, 5, 10, 15],
                  regressor__n_estimators=[10, 50, 100, 150, 200],
                  regressor__max_features=["auto", "log2", "sqrt"])
                  #regressor__C=[0.001, 0.01, 0.1, 1, 10],
                  #regressor__gamma=[0.001, 0.01, 0.1, 1])
                  #regressor__class_weight=[{0:1, 1:2},{0:1, 1:3}, {0:1, 1:4}, {0:1, 1:5}])

#x = data.iloc[:, :-1]
#y = data.iloc[:, -1]

#X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.33, random_state=4, shuffle=True)



grid_search = GridSearchCV(pipeline, param_grid=param_grid, verbose=10, n_jobs=4, cv=5, refit=True)
grid_search.fit(X_train, y_train)
estimator = grid_search.best_estimator_
#predictions = estimator.predict(X_test)
#probability = estimator.predict_proba(X_test)[:,1]

#y_pred = pd.DataFrame(predictions, columns=['Predictions']).to_csv('predictions.csv', index=False, header=False)
#prob = pd.DataFrame(probability, columns=['Probability']).to_csv('prob_ACRinc_di_SVC.csv', index=False, header=False)
#y_true = pd.DataFrame(y_test).to_csv('True.csv', index=False, header=False)

with open('model.pkl','wb') as f:
    pickle.dump(estimator,f)

# load
#with open('model.pkl', 'rb') as f:
#    clf2 = pickle.load(f)
#
#clf2.predict(X[0:1])

#print("Coef")
#print(grid_search.best_estimator_.named_steps["regressor"].coef_)
#print("Ranking")
#print(grid_search.best_estimator_.named_steps["RFE"].ranking_)
