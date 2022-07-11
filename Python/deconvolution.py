import pandas as pd
from pycaret.datasets import get_data
from pycaret.classification import *
import numpy as np
from pycaret.utils import check_metric
from sklearn.metrics import classification_report, confusion_matrix

home_directory = "/home/ottoraik"
target_vector = pd.read_csv(home_directory + "/Dropbox/testproject/Datasets/Target_vector_binary.tsv",sep  = "\t", index_col=0)
dataset = pd.read_csv(home_directory + "/Dropbox/testproject/Datasets/Deconvolution_RepSet_S84.tsv",sep  = "\t")
#dataset = pd.read_csv(home_directory + "/Dropbox/testproject/Datasets/Deconvolution_RepSet_S84.946genes.training_set.tsv",sep  = "\t")

dataset.columns
#dataset = dataset.reset_index(drop=True,inplace=True)
dataset.index
target_vector.index
target_vector.reset_index(drop=True,inplace=True)
target_vector.index

#dataset = pd.DataFrame(dataset.transpose())
dataset = pd.concat( [dataset ,target_vector],axis = 1)
#dataset = dataset.rename(columns = {"truth_vec_binary": "Grading_binary"})
dataset.columns

dataset = dataset.drop(columns =["truth_vec_binary","Sample","Grading","Grading_binary_numeric"])
dataset.columns

'indicator_vec  = dataset[["Grading"]].isin(["G1","G2"])'
'dataset.loc[5,"Grading"]'
'dataset[["Grading"]] = ["G1_G2" if indicator_vec == True]'

'dataset[["Grading"]].query("Grading  in ["G1","G2"]")'

#selector = ["alpha","beta","ductal","P_value","Correlation","RMSE","Grading_binary"]

exp_clf102 = setup(
    data = dataset,
    target = 'Grading_binary',
    normalize = True,
    transformation = True,
    ignore_low_variance = True,
    remove_multicollinearity = False,
    numeric_imputation = "median",
    log_experiment = True, experiment_name = 'Deco_1',
    feature_selection = True)

best_model = compare_models(sort='Accuracy')

del rf
rf = create_model('rf',fold=10)
result_rf = pull()

tuned_rf = tune_model(rf, optimize = 'Accuracy',fold=10)
result_rf = pull()
#plot_model(tuned_rf, plot = 'auc')
#plot_model(tuned_rf, plot = 'rp')
#plot_model(tuned_rf, plot='feature')
#plot_model(tuned_rf, plot = 'confusion_matrix')

final_rf = finalize_model(tuned_rf)
result_final_model = pull()

validation_predictions = predict_model(final_rf, data=training_data)
#check_metric(unseen_predictions['Grading_binary'], dataset['Grading_binary'], metric = 'Accuracy')
confusion_matrix(training_data['Grading_binary'],validation_predictions['Label'])
results = classification_report(training_data['Grading_binary'],validation_predictions['Label'])

import autosklearn.classification
cls = autosklearn.classification.AutoSklearnClassifier()
cls.fit(X_train, y_train)
predictions = cls.predict(X_test)
