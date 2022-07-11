import pandas as pd
import pickle
from pycaret.datasets import get_data
from pycaret.classification import *
import numpy as np
from pycaret.utils import check_metric
from sklearn.metrics import classification_report, confusion_matrix

#target_vector = pd.read_csv("/home/ottoraik/Dropbox/testproject/Target_vector_drug_response.tsv",sep  = "\t", index_col=0)
#target_vector = target_vector[["V2"]]
target_vector = pd.read_csv("/home/ottoraik/Dropbox/testproject/IC50.tsv",sep  = "\t", index_col=0)
#dataset = pd.read_csv("/home/ottoraik/Dropbox/testproject/Data_9461.Counts.HGNC.tsv",sep  = "\t", index_col = 0)
dataset = pd.read_csv("/home/ottoraik/Dropbox/testproject/IC50_exp.tsv",sep  = "\t", index_col = 0)

#dataset = pd.DataFrame(dataset.transpose())
target_vector = target_vector.set_index(dataset.index)

dataset_train = pd.concat( [dataset,target_vector],axis = 1)
#dataset = dataset.reset_index(drop=True,inplace=True)
#dataset_train = dataset_train.rename(columns = {"V2": "Drug_treatment"})

#dataset = dataset.drop(columns =["Grading"])
#dataset = dataset.drop(columns =["Grading_binary"])
#dataset = dataset.drop(columns =["Grading_binary_numeric"])
#dataset = dataset.drop(columns =["Sample"])

training_data = dataset_train
#training_data   = dataset.sample(frac=0.80, random_state = 123)
#validation_data = dataset.drop(training_data.index)
#dataset = dataset.reset_index(drop=True,inplace=True)

'indicator_vec  = dataset[["Grading"]].isin(["G1","G2"])'
'dataset.loc[5,"Grading"]'
'dataset[["Grading"]] = ["G1_G2" if indicator_vec == True]'

'dataset[["Grading"]].query("Grading  in ["G1","G2"]")'

#selector = ["alpha","beta","ductal","P_value","Correlation","RMSE","Grading_binary"]

exp_clf102 = setup(
    data = training_data,
    target = ['IC50'],
    normalize = True,
    transformation = True,
    ignore_low_variance = False,
    remove_multicollinearity = True,
    numeric_imputation = "median",
    log_experiment = True,
    experiment_name = 'IC50',
    feature_selection = True)

best_model = compare_models(sort='Accuracy')

del rf
rf = create_model('svm',fold=10)
result_rf = pull()

tuned_rf = tune_model(rf, optimize = 'Accuracy',fold=10)
result_rf = pull()
#plot_model(tuned_rf, plot = 'auc')
#plot_model(tuned_rf, plot = 'rp')
plot_model(final_rf, plot='feature')
#plot_model(tuned_rf, plot = 'confusion_matrix')

final_rf = finalize_model(tuned_rf)
result_final_model = pull()

#pickle.dump(final_rf, open("/home/ottoraik/Dropbox/testproject/svm_linear.sav", 'wb'))
#pickleFile = open("/home/ottoraik/Dropbox/testproject/svm_linear.sav", 'rb')
final_rf = pickle.load(pickleFile)

validation_predictions = predict_model(final_rf, data=training_data)
#check_metric(unseen_predictions['Grading_binary'], dataset['Grading_binary'], metric = 'Accuracy')
confusion_matrix(training_data['Drug_treatment'],validation_predictions['Label'])
results = classification_report(training_data['Drug_treatment'],validation_predictions['Label'])

features = pd.DataFrame({'Feature': get_config('X_train').columns, 'Value' : abs(final_rf.coef_[0])}).sort_values(by='Value', ascending=False)
pd.DataFrame.to_csv(features, "/home/ottoraik/Dropbox/testproject/features.csv")

validation_dataset = pd.read_csv("/home/ottoraik/Dropbox/testproject/Schmitz_S2546.HGNC.tsv",sep  = "\t")
#validation_dataset = validation_dataset.set_index('Unnamed: 0')
#validation_dataset = validation_dataset.transpose()

validation_predictions = predict_model(final_rf, data=validation_dataset)

pd.DataFrame.to_csv(validation_predictions[['Label']], "/home/ottoraik/Dropbox/testproject/Schmit.csv", sep ="\t")