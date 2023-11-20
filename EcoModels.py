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
parser.add_argument("--rf_trees", help="Number of trees to grow in each RF", default=100, type = int)
parser.add_argument("--min_pred_prev", help="Minimum number of non-zero values required for variables in the preditor DF to be include in models", default=1, type = int)
parser.add_argument("--min_resp_prev", help="Minimum number of non-zero values required for variables in the response DF to be include in models", default=1, type = int)
parser.add_argument("--k_cv", help="Value of k t be used for k-fold cross validations. Skipped if k = 0", default=0, type = int)
args = parser.parse_args()


def central():
    #Format the datasets fto be used by the models
    (full_df,pred_df,resp_df,valid_pred_var_names,valid_resp_var_names) = format_datasets()
    #Build ANN models
    if (args.build_ann == True):
        #Format the datasets for ANN models (i.e. any columns with NA values in any samples are removed)
        full_MOR = build_ann(pred_df=pred_df,resp_df=resp_df)
        full_perfo_info_df = calc_performance_metrics(mor_model=full_MOR, pred_df=pred_df,resp_df=resp_df, print_preds=True)
        full_perfo_info_df["Dataset"] = "Full"
        output_dataframe_file = args.prefix + "_Full_Models_Performance_Info.tsv"
        full_perfo_info_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    #Build RF models
    if (args.build_rf == True):
        rf_full_df = full_df
        #Run k-fold cross validation if specified
        if (args.k_cv > 0):
            print(f"Running {args.k_cv}-fold cross validation")
            kf = KFold(n_splits=args.k_cv,shuffle=True,random_state=42)
            fold_idx = 0
            perfo_dfs = []
            for train_index, test_index in kf.split(rf_full_df):
                fold_idx = fold_idx + 1
                print(f"Running fold {fold_idx}")
                #Split into training and test sets
                full_cv_train_df = rf_full_df.iloc[train_index]
                full_cv_test_df = rf_full_df.iloc[test_index]
                #Split into predictor and response variables
                cv_train_pred_df = full_cv_train_df[valid_pred_var_names]
                cv_train_resp_df = full_cv_train_df[valid_resp_var_names]
                cv_MOR = build_rf(pred_df=cv_train_pred_df,resp_df=cv_train_resp_df,rf_trees=args.rf_trees)
                #Calculate performance metrics of the CV models on Training set
                cv_train_perfo_info_df = calc_performance_metrics(mor_model=cv_MOR, pred_df=cv_train_pred_df,resp_df=cv_train_resp_df)
                cv_train_perfo_info_df["Fold"] = fold_idx
                cv_train_perfo_info_df["Dataset"] = "Training"
                #Calculate performance metrics of the CV models on Test set
                cv_test_pred_df = full_cv_test_df[valid_pred_var_names]
                cv_test_resp_df = full_cv_test_df[valid_resp_var_names]
                cv_test_perfo_info_df = calc_performance_metrics(mor_model=cv_MOR, pred_df=cv_test_pred_df,resp_df=cv_test_resp_df)
                cv_test_perfo_info_df["Fold"] = fold_idx
                cv_test_perfo_info_df["Dataset"] = "Test"
                #Merge the two dataframes
                cv_perfo_info_df = pd.concat([cv_train_perfo_info_df, cv_test_perfo_info_df], axis=0)
                perfo_dfs.append(cv_perfo_info_df)
            #Merge the performance dataframes from each fold
            full_cv_perfo_info_df = pd.concat(perfo_dfs, axis=0)
            output_dataframe_file = args.prefix + "_CV_Models_Performance_Info.tsv"
            full_cv_perfo_info_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
        #Build a model using the full dataset
        full_MOR = build_rf(pred_df=rf_full_df[valid_pred_var_names],resp_df=rf_full_df[valid_resp_var_names],rf_trees=args.rf_trees)
        full_perfo_info_df = calc_performance_metrics(mor_model=full_MOR, pred_df=rf_full_df[valid_pred_var_names],resp_df=rf_full_df[valid_resp_var_names], print_preds=True)
        full_perfo_info_df["Dataset"] = "Full"
        output_dataframe_file = args.prefix + "_Full_Models_Performance_Info.tsv"
        full_perfo_info_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')

def format_datasets():
    if args.resp_pred_table:
        #Warn user if the resp_pred_file is provide dby the predictor and response variables are empty, then die
        if (args.predictor_vars == None) or (args.response_vars == None):
            print("WARNING: When porividing the resp_pred_file both predictor AND response variables must be specified")
            exit()
        #For simplicity, the predictor and response variables are loaded and processed separately then merged
        (resp_df, valid_resp_var_names) = prepare_single_df(type="response",df_file=args.resp_pred_table,min_var_prev=args.min_resp_prev,idx_var=args.response_index,transpose=args.transpose_response,z_transform=args.z_transform,valid_var_names=args.response_vars)
        (pred_df, valid_pred_var_names) = prepare_single_df(type="predictor",df_file=args.resp_pred_table,min_var_prev=args.min_pred_prev,idx_var=args.predictor_index,transpose=args.transpose_predictors,z_transform=args.z_transform,valid_var_names=args.predictor_vars)
    else:
        (resp_df, valid_resp_var_names) = prepare_single_df(type="response",df_file=args.response_table,min_var_prev=args.min_resp_prev,idx_var=args.response_index,transpose=args.transpose_response,z_transform=args.z_transform,valid_var_names=args.response_vars)
        (pred_df, valid_pred_var_names) = prepare_single_df(type="predictor",df_file=args.predictor_table,min_var_prev=args.min_pred_prev,idx_var=args.predictor_index,transpose=args.transpose_predictors,z_transform=args.z_transform,valid_var_names=args.predictor_vars)
    
    full_df = pd.concat([resp_df, pred_df], axis=1, join='inner')
    print(f"Shape of full DF: {full_df.shape}")
    output_dataframe_file = args.prefix + "_Full_Processed_Data_Frame.tsv"
    full_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    pred_df = full_df[valid_pred_var_names]
    resp_df = full_df[valid_resp_var_names]
    return(full_df,pred_df,resp_df,valid_pred_var_names,valid_resp_var_names)

@ignore_warnings(category=ConvergenceWarning)
def build_ann(resp_df=None, pred_df=None):
    #Perform hyper parameter tuning using grid search
    param_grid={'estimator__learning_rate_init': [0.001, 0.005, 0.01], 'estimator__max_iter': [100], 'estimator__hidden_layer_sizes': [2,3,5]}
    hpt_ann = MLPRegressor(random_state=42)
    hpt_MOR = MultiOutputRegressor(hpt_ann,n_jobs=1)
    print(f"Performing ANN hyperparameter tuning")
    gsMOR = GridSearchCV(estimator=hpt_MOR, param_grid=param_grid, refit=True, cv=KFold(n_splits=5), verbose=3, n_jobs=args.threads, return_train_score=True, scoring="r2") #cv=5, scoring="r2"neg_mean_squared_error
    gsMOR.fit(pred_df, resp_df)
    print(f"Best ANN hyperparameters: {gsMOR.best_params_}")
    
    hpt_results_df = pd.DataFrame(gsMOR.cv_results_)
    output_dataframe_file = args.prefix + "_ANN_Hyper_Parameter_Tuning_Metrics.tsv"
    hpt_results_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    #Buiild a model for each response variable using the MultiOutputRegressor function
    #print(f"Training ANN models")
    #Set up ANN. Hyperparams are the same for all the output variables
    #ann = MLPRegressor(hidden_layer_sizes=10, max_iter=1000, learning_rate_init=0.01, random_state=42)
    # define the direct multioutput MOR model
    #MOR = MultiOutputRegressor(ann,n_jobs=args.threads)
    #MOR.fit(pred_df, resp_df)
    #Calculate performance metrics of the models
    best_MOR = gsMOR.best_estimator_
    #Print the number of hidden layers in the first object in best_MOR
    print(f"Total Number of layers (input + hidden + output) in best Model: {best_MOR.estimators_[0].n_layers_}")
    return(best_MOR)

def calc_performance_metrics(mor_model=None, pred_df=None, resp_df=None, print_preds=False):
    print("Calculating model performance metrics")
    model_info = defaultdict(dict)
    predictions = mor_model.predict(pred_df)
    resp_var_names = resp_df.columns
    predictions_df = pd.DataFrame(predictions, columns=resp_var_names, index=resp_df.index)
    predictions_df.index.name = resp_df.index.name
    #predictions_df["Sample"] = pred_df.index
    # print(f"Summary of model predictions:")
    # print(predictions_df.describe())
    #predictions_df.index.name = 'Response_Variable'
    #print(predictions_df.describe())
    for rvar in resp_var_names:
        model_info[rvar]["PearsonR"] = scipy.stats.pearsonr(predictions_df[rvar], resp_df[rvar])[0]
        model_info[rvar]["RMSE"] = np.sqrt(np.mean((predictions_df[rvar] - resp_df[rvar].values.ravel())**2))
        #print(f"Response Variable: {rvar}, PearsonR: {model_info[rvar]['PearsonR']}, RMSE: {model_info[rvar]['RMSE']}")
    model_info_df = pd.DataFrame.from_dict(model_info)
    model_info_df = model_info_df.transpose()
    model_info_df.index.name = 'Response_Variable'

    if print_preds == True:
        #Merge the predictions with the true values of the response variables in resp_df indicating which is each
        predictions_df["Data_Type"] = "Predicted"
        resp_df["Data_Type"] = "Measured"
        predictions_df = pd.concat([predictions_df, resp_df], axis=0)
        output_dataframe_file = args.prefix + "_Full_Models_Predicted_vs_Measured_Values.tsv"
        predictions_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    
    return(model_info_df)

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
    return(MOR)


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
