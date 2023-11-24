#! /usr/bin/env python3
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from collections import defaultdict
import scipy.stats
import pickle
#from scipy.stats import zscore
import pandas as pd
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--resp_pred_table", help=".tsv format file containing the both the predictor and response variables", type = str)
parser.add_argument("--response_table", help=".tsv format file containing the response variables", type = str)
parser.add_argument("--predictor_table", help=".tsv format file containing the predictor variables", type = str)
parser.add_argument("--predictor_vars", help="list of variables in the predictor table to be used when building models. If none are specified, all variables in the table will be used", nargs="*", type = str)
parser.add_argument("--response_vars", help="list of variables in the response table to be used when building models. If none are specified, all variables in the table will be used", nargs="*", type = str)
parser.add_argument("--transpose_predictors", help="Flag to transpose the predictor_table file", default=False, type = bool)
parser.add_argument("--transpose_response", help="Flag to transpose the response_table file", default=False, type = bool)
parser.add_argument("--transpose_resp_pred", help="Flag to transpose the response_table file", default=False, type = bool)
parser.add_argument("--response_index", help="Response variables tables index column name. If not provided assumed to be the first column", type = str)
parser.add_argument("--predictor_index", help="Predictor variables tables index column name. If not provided assumed to be the first column", type = str)
parser.add_argument("--resp_pred_index", help="Predictor variables tables index column name. If not provided assumed to be the first column", type = str)
parser.add_argument("--sep_table", help="Input tables separator charcter", nargs=1, default="\t", type = str)
parser.add_argument("--prefix", help="Prefix of output file names", default="EcoModels", type = str)
parser.add_argument("--z_transform", help="Flag to Z-transform values of predictor and response variables before training models", default=False, type = bool)
parser.add_argument("--threads", help="Number of threads to be used during analysis", default=1, type = int)
parser.add_argument("--build_ann", help="Flag to run the ANN module", default=False, type = bool)
parser.add_argument("--build_rf", help="Flag to run the Random Forest module", default=False, type = bool)
parser.add_argument("--min_pred_prev", help="Minimum number of non-zero values required for variables in the preditor DF to be include in models", default=1, type = int)
parser.add_argument("--min_resp_prev", help="Minimum number of non-zero values required for variables in the response DF to be include in models", default=1, type = int)
parser.add_argument("--run_hpt", help="Flag to run hyperparameter tuning with grid search cross validation to find best HP to be used for training a final model", default=False, type = bool)
parser.add_argument("--hpt_rf_trees", help="Number of trees to grow in each RF", default=[100,500,1000], type = int, nargs="+")
parser.add_argument("--hpt_k_cv", help="Value of k t be used for k-fold cross validations durign hyperparameter tuning. Skipped if k = 0", default=5, type = int)
parser.add_argument("--hpt_ann_n_neurons", help="Values of number of neurons in the hidden layer to be evaluated during ANN hyperparameter tuning", default=[2,3,5,10], type = int,  nargs="+")
parser.add_argument("--hpt_ann_learning_rate_init", help="Values of number of initial learning rate to be evaluated during ANN hyperparameter tuning", default=[0.001, 0.005, 0.01, 0.05], type = float,  nargs="+")
parser.add_argument("--hpt_ann_max_iter", help="Values of maximum number of iterations to be evaluated during ANN hyperparameter tuning", default=[50,100,500], type = int,  nargs="+")
args = parser.parse_args()


def central():
    #Format the datasets to be used by the models any columns with NA values in any samples are removed)
    (full_df,pred_df,resp_df,valid_pred_var_names,valid_resp_var_names) = format_datasets()
    #Build ANN models
    if (args.build_ann == True):
        full_MOR = build_model(pred_df=pred_df,resp_df=resp_df,model_type="ANN")
        #full_MOR = build_ann(pred_df=pred_df,resp_df=resp_df)
        (performance_df, predictions_df) = calc_performance_metrics(mor_model=full_MOR, pred_df=pred_df,resp_df=resp_df, model_type="ANN", dataset="Full", print_perfo=True, print_preds=True)
    #Build RF models
    if (args.build_rf == True):
        #Build a model using the full dataset
        full_MOR = build_model(pred_df=pred_df,resp_df=resp_df,model_type="RF")
        (performance_df, predictions_df) = calc_performance_metrics(mor_model=full_MOR, pred_df=pred_df,resp_df=resp_df, model_type="RF", dataset="Full", print_perfo=True, print_preds=True)

def calc_performance_metrics(mor_model=None, pred_df=None, resp_df=None, model_type=None, dataset="NA", print_perfo=False, print_preds=False):
    print("Calculating model performance metrics")
    model_info = defaultdict(dict)
    #Generarte predictions for the response variables based on the output of the received model
    predictions = mor_model.predict(pred_df)
    resp_var_names = resp_df.columns
    predictions_df = pd.DataFrame(predictions, columns=resp_var_names, index=resp_df.index)
    predictions_df.index.name = resp_df.index.name
    #Calculate model performance metrics for each response variable in the MOR model
    for rvar in resp_var_names:
        model_info[rvar]["PearsonR"] = scipy.stats.pearsonr(predictions_df[rvar], resp_df[rvar])[0]
        model_info[rvar]["RMSE"] = np.sqrt(np.mean((predictions_df[rvar] - resp_df[rvar].values.ravel())**2))
        #print(f"Response Variable: {rvar}, PearsonR: {model_info[rvar]['PearsonR']}, RMSE: {model_info[rvar]['RMSE']}")
    performance_df = pd.DataFrame.from_dict(model_info)
    performance_df = performance_df.transpose()
    performance_df.index.name = 'Response_Variable'
    if (print_perfo == True):
        performance_df["Dataset"] = dataset
        output_dataframe_file = args.prefix + model_type + dataset + "_Models_Performance_Info.tsv"
        performance_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    
    #Merge the predictions with the expected/true values of the response variables in resp_df indicating which is each. This must be done after the performance metrics are calculated so that the DF maintain the correc dimensions
    predictions_df["Data_Type"] = "Predicted"
    resp_df["Data_Type"] = "Measured"
    predictions_df = pd.concat([predictions_df, resp_df], axis=0)
    if (print_preds == True):
        output_dataframe_file = args.prefix + model_type + dataset  + "_Models_Predicted_vs_Measured_Values.tsv"
        predictions_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    
    return(performance_df,predictions_df)

def build_model(resp_df=None, pred_df=None, model_type=None):
    #Instantiate a model deppending on the type of model to be built
    if (model_type == "ANN"):
        #Set up an ANN model with default hyperparamters= {'learning_rate_init': 0.001, 'max_iter': 100, 'hidden_layer_sizes': (2,)}
        model = MLPRegressor(random_state=42,learning_rate_init=0.01, max_iter=100, hidden_layer_sizes=(5,))
        param_grid={'estimator__learning_rate_init': args.hpt_ann_learning_rate_init, 'estimator__max_iter': args.hpt_ann_max_iter, 'estimator__hidden_layer_sizes': args.hpt_ann_n_neurons}
    elif (model_type == "RF"):
        model = RandomForestRegressor(random_state=42, n_jobs=1, n_estimators=100)
        param_grid={'estimator__n_estimators': args.hpt_rf_trees}
    #Instantiate a MultiOutputRegressor object to build a model for each response variable
    MOR = MultiOutputRegressor(model,n_jobs=1)
    #Instantiate a GridSearchCV object to perform hyperparameter tuning
    if (args.run_hpt == True):
        #Set combinations of hyperparameters to be tested
        HPT = GridSearchCV(estimator=MOR, param_grid=param_grid, refit=True, cv=KFold(n_splits=args.hpt_k_cv, shuffle=True), verbose=0, n_jobs=args.threads, return_train_score=True, scoring="neg_mean_squared_error")
        #Perform the grid search, obtain best parameters, and fit the final model
        HPT.fit(pred_df, resp_df)
        hpt_results_df = pd.DataFrame(HPT.cv_results_)
        output_dataframe_file = args.prefix + model_type +"_Hyper_Parameter_Tuning_Metrics.tsv"
        hpt_results_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
        best_MOR = HPT.best_estimator_
    else:
        #Fit the single model
        MOR.fit(pred_df, resp_df)
        best_MOR = MOR
    #Save model to file with pickle.dump
    output_model_file = args.prefix + "_" + model_type +"_MOR_Model.pkl"
    pickle.dump(best_MOR, open(output_model_file, 'wb'))
    return(best_MOR)

def build_ann(resp_df=None, pred_df=None):
    #Instantiate ANN
    ann = MLPRegressor(random_state=42)
    #Buiild a model for each response variable using the MultiOutputRegressor function
    MOR = MultiOutputRegressor(ann,n_jobs=1)
    print(f"Performing ANN hyperparameter tuning")
    #Set combinations of hyperparameters to be tested
    param_grid={'estimator__learning_rate_init': args.hpt_ann_learning_rate_init, 'estimator__max_iter': args.hpt_ann_max_iter, 'estimator__hidden_layer_sizes': args.hpt_ann_n_neurons}
    #Set up the grid search
    HPT = GridSearchCV(estimator=MOR, param_grid=param_grid, refit=True, cv=KFold(n_splits=args.hpt_k_cv, shuffle=True), verbose=0, n_jobs=args.threads, return_train_score=True, scoring="neg_mean_squared_error") #cv=5, scoring="r2"neg_mean_squared_error
    #Perform the grid search, obtain best parameters, and fit the final model
    HPT.fit(pred_df, resp_df)
    #print(f"Best ANN hyperparameters: {HPT.best_params_}")
    #The best score is actually defined as the mean_test_score over all folds. This is the score that is used to select the best parameters before refitting a final
    #print(f"Best Score: {HPT.best_score_}")
    #Builda DF from the results of the grid search and print it to a file
    hpt_results_df = pd.DataFrame(HPT.cv_results_)
    output_dataframe_file = args.prefix + "_ANN_Hyper_Parameter_Tuning_Metrics.tsv"
    hpt_results_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')

    #Obtain the best models from the grid search after refitting
    best_MOR = HPT.best_estimator_
    #Print the number of hidden layers in the first object in best_MOR
    #print(f"Total Number of layers (input + hidden + output) in best Model: {best_MOR.estimators_[0].n_layers_}")
    #Save model to file with pickle.dump
    output_model_file = args.prefix + "_Best_ANN_MOR_Model.pkl"
    pickle.dump(best_MOR, open(output_model_file, 'wb'))

    return(best_MOR)

def build_rf(pred_df=None, resp_df=None, rf_trees=100):
    """Build a Random Forest model for each response variable according to the user specifications"""
    #Buiild a model for each response variable.
    MOR = MultiOutputRegressor(RandomForestRegressor(n_estimators=rf_trees,random_state=42,n_jobs=args.threads),n_jobs=1)
    MOR.fit(pred_df, resp_df)
    #rf.fit(pred_df, resp_df)
    #Iterate over each generated RF model
    #Could also do it with MOR.feature_names_in_
    resp_var_names = resp_df.columns
    idx = 0
    imp_dfs_list = []
    for rf in MOR.estimators_:
        #Calculate the relative importance of predictors
        rname = resp_var_names[idx]
        #Calculate relative importance of predictor variables
        pred_imps = rf.feature_importances_
        model_imp_df = pd.DataFrame(pred_imps)
        model_imp_df.index.name = 'Predictor_Index'
        model_imp_df.rename(columns={0: rname},inplace=True)
        model_imp_df['Predictor_Variable'] = pred_df.columns
        model_imp_df.set_index('Predictor_Variable',inplace=True)
        imp_dfs_list.append(model_imp_df)
        #output_dataframe_file = args.prefix +f"_Response_Variable_{rname}"+ "_Predictor_Importance.tsv"
        #model_imp_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
        idx = idx + 1
    full_imp_df = pd.concat(imp_dfs_list,axis=1)
    output_dataframe_file = args.prefix + "_Predictor_Importance.tsv"
    #full_imp_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    #Save model to file with pickle.dump
    output_model_file = args.prefix + "_Best_RF_MOR_Model.pkl"
    pickle.dump(MOR, open(output_model_file, 'wb'))
    return(MOR)
    
def format_datasets():
    if args.resp_pred_table:
        #Warn user if the resp_pred_file is provide dby the predictor and response variables are empty, then die
        if (args.predictor_vars == None) or (args.response_vars == None):
            print("WARNING: When porividing the resp_pred_file both predictor AND response variables must be specified")
            exit()
        #For simplicity, the predictor and response variables are loaded and processed separately then merged
        resp_file = args.resp_pred_table
        pred_file = args.resp_pred_table
    else:
        resp_file = args.response_table
        pred_file = args.predictor_table
    #Load and process the predictor and response variables
    (resp_df, valid_resp_var_names) = prepare_single_df(type="response",df_file=resp_file,min_var_prev=args.min_resp_prev,idx_var=args.response_index,transpose=args.transpose_response,z_transform=args.z_transform,valid_var_names=args.response_vars)
    (pred_df, valid_pred_var_names) = prepare_single_df(type="predictor",df_file=pred_file,min_var_prev=args.min_pred_prev,idx_var=args.predictor_index,transpose=args.transpose_predictors,z_transform=args.z_transform,valid_var_names=args.predictor_vars)
    
    full_df = pd.concat([resp_df, pred_df], axis=1, join='inner')
    print(f"Shape of full DF: {full_df.shape}")
    output_dataframe_file = args.prefix + "_Full_Processed_Data_Frame.tsv"
    full_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    pred_df = full_df[valid_pred_var_names]
    resp_df = full_df[valid_resp_var_names]
    return(full_df,pred_df,resp_df,valid_pred_var_names,valid_resp_var_names)








def prepare_single_df(type=None,df_file=None,min_var_prev=1,idx_var=None,transpose=False,z_transform=False,valid_var_names=None):
    print(f"Reading in {type} data from: {df_file}")
    if idx_var:
        print(f"Using: ",idx_var," as index column")
    else:
        idx_var = 0
        print(f"Using first column as index")
    full_df = pd.read_csv(df_file,sep=args.sep_table,index_col=idx_var,header=0)
    print(f"Shape of full DF: {full_df.shape}")
    #Transpose DF if specified
    if transpose == True:
        print("Transposing data frame")
        full_df = full_df.transpose()
        print(f"Shape of transposed DF: {full_df.shape}")
    #If no response variables are specified, keep all variables in the table
    all_var_names = full_df.columns
    #If a specific list of variables is provided, keep only those
    sub_df = full_df
    if valid_var_names:
        print(f"Keeping only specified variables: {valid_var_names}")
        sub_df = full_df[valid_var_names]
        print(f"Shape of DF after removing non-selected variables: {sub_df.shape}")  
    else:
        #Keep only the user specified variables
        print("Using all variables")
        valid_var_names = all_var_names  
    #remove any rows with missing values from sub_df
    print("Removing rows with missing values")
    sub_df.dropna(inplace=True, axis=0, how='any')
    print(f"Shape of DF after removing missing values: {sub_df.shape}")
    #Count prevalence of non-zero values for variable column in the DF
    var_prev = sub_df.astype(bool).sum(axis=0)
    #Keep only the variables for which the required minimum prevalence is not met
    valid_var_names = var_prev[var_prev >= min_var_prev].index
    sub_df = sub_df[valid_var_names]
    print(f"{len(valid_var_names)} variables kept after prevalence filtering")
    #Calculated standard deviation of each variable in sub_df
    var_stdev = sub_df.std(axis=0)
    #Remove variables with zero standard deviation
    valid_var_names = var_stdev[var_stdev > 0].index
    sub_df = sub_df[valid_var_names]
    print(f"{len(valid_var_names)} variables kept after standard deviation filtering")
    #Z-transfor the DF if specified
    if (z_transform == True):
        print("Z-transforming data")
        sub_df = sub_df.apply(scipy.stats.zscore)
        print(f"Shape of DF after Z-transforming values: {sub_df.shape}")
    #Return the formatted DF and the list with the names of valid variables
    valid_var_names = sub_df.columns
    print(f"Summary of DF after pre-processing:")
    print(sub_df.describe())
    return(sub_df,valid_var_names)


     
    

   




central()
