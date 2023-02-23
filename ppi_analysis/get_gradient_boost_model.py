#!/usr/bin/python

from numpy import arange
from pandas import read_csv
from sklearn.model_selection import RepeatedKFold
from pandas import read_csv
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import plot_roc_curve
import random
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu

from matplotlib import pyplot

####### Gradient Boosting Classifier ####
from sklearn.ensemble import GradientBoostingClassifier


import time
start_time = time.time()

model = GradientBoostingClassifier(n_estimators=50, learning_rate=0.4, max_depth=1, random_state=0)

# Load train dataset
train_url =  "/data/timonaj/cancer_as_wound/ppi_analysis/machine_learning/STRING_S_EXP_total_fraction_total.csv"
dataframe_train = read_csv(train_url, header=0)
data_train = dataframe_train.values
X_train, Y_train = data_train[:, 1:-1], data_train[:, -1]


model.fit(X_train, Y_train)

# save the model to disk

filename = 'finalized_STRING_model.sav'

pickle.dump(model, open(filename, 'wb'))

print("--- %s seconds ---" % (time.time() - start_time))