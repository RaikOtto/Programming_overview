import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from sklearn.metrics import confusion_matrix, classification_report
from pycaret.datasets import get_data
from itertools import product
from pycaret.regression import *
from pycaret.utils import check_metric

target_vector = pd.read_csv("/home/ottoraik/Dropbox/testproject/IC50_log.csv",sep  = "\t", index_col=0)
dataset = pd.read_csv("/home/ottoraik/Dropbox/testproject/IC50_exp.csv",sep  = ",", index_col = 0)

target_vector = target_vector.set_index(dataset.index)

dataset_train = pd.concat( [dataset,target_vector],axis = 1)
training_data = dataset_train

exp_clf102 = setup(
    data = training_data,
    target = 'IC50',
    normalize = True,
    transformation = True,
    ignore_low_variance = False,
    remove_multicollinearity = True,
    numeric_imputation = "median",
    log_experiment = True,
    experiment_name = 'IC50',
    feature_selection = True)

best_model = compare_models(sort='RMSE')
result_rf = pull()
result_rf

lr = create_model('ada')
plot_model(lr, plot='feature')
plot_model(lr)

predict_model(lr);
results = pull()
results

final_lr = finalize_model(lr)
importances = final_lr.feature_importances_

pickle.dump(final_lr, open("/home/ottoraik/Dropbox/testproject/ada_ic50_log_regression.sav", 'wb'))

validation_predictions = predict_model(final_lr, data=training_data)
predictions = validation_predictions[["Label"]]
original = training_data[["IC50"]]
pd.concat( [predictions,original],axis = 1)

#features = pd.DataFrame({'Feature': get_config('X_train').columns, 'Value' : abs(final_lr.coef_[0])}).sort_values(by='Value', ascending=False)
features = pd.DataFrame({'Feature': get_config('X_train').columns, 'Value' : abs(final_lr.feature_importances_)}).sort_values(by='Value', ascending=False)
pd.DataFrame.to_csv(get_config('X_train').columns, "/home/ottoraik/Dropbox/testproject/IC50_log_regression.features.csv",sep ="\t")

validation_dataset = pd.read_csv("/home/ottoraik/Dropbox/testproject/chapuy_ic50_exp.csv",sep  = "\t")
validation_dataset = pd.DataFrame(validation_dataset.transpose())
validation_dataset
validation_dataset = validation_dataset.set_index(validation_dataset.index)
validation_dataset.index
validation_predictions = predict_model(final_lr, data=validation_dataset)
pd.DataFrame.to_csv(validation_predictions[['Label']], "/home/ottoraik/Dropbox/testproject/Chapuy_IC50.csv", sep ="\t")