#! /usr/bin/env python3
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.ensemble import RandomForestRegressor
from collections import defaultdict
import scipy.stats
import pandas as pd
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--response_table", help=".tsv format file containing the response variables", required=True, type = str)
parser.add_argument("--predictor_table", help=".tsv format file containing the predictor variables", required=True, type = str)
parser.add_argument("--predictor_vars", help="list of variables in the predictor table to be used when building models. If none are specified, all variables in the table will be used", nargs="*", type = str)
parser.add_argument("--response_vars", help="list of variables in the response table to be used when building models. If none are specified, all variables in the table will be used", nargs="*", type = str)
parser.add_argument("--transpose_predictors", help="Flag to transpose the predictor_table file", default=False, type = bool)
parser.add_argument("--transpose_response", help="Flag to transpose the response_table file", default=False, type = bool)
parser.add_argument("--response_index", help="Response variables tables index column name. If not provided assumed to be the first column", type = str)
parser.add_argument("--predictor_index", help="Predictor variables tables index column name. If not provided assumed to be the first column", type = str)
parser.add_argument("--sep_table", help="Input tables separator charcter", nargs=1, default="\t", type = str)
parser.add_argument("--prefix", help="Prefix of output file names", default="EcoModels", type = str)
parser.add_argument("--threads", help="Number of threads to be used during analysis", default=1, type = int)
parser.add_argument("--build_ann", help="Flag to run the ANN module", default=False, type = bool)
parser.add_argument("--build_rf", help="Flag to run the Random Forest module", default=False, type = bool)
args = parser.parse_args()


def central():
    #Format the datasets
    (resp_df, pred_df) = format_datasets()
    #Build ANN models
    if (args.build_ann == True):
        build_ann(resp_df=resp_df, pred_df=pred_df)
    #Build RF models
    if (args.build_rf == True):
        build_rf(resp_df=resp_df, pred_df=pred_df)

def build_rf(resp_df=None, pred_df=None):
    #Buiild a model for each response variable.
    rf = RandomForestRegressor(n_estimators=1000,random_state=42,n_jobs=args.threads)
    rf.fit(pred_df, resp_df)
    pred_imps = rf.feature_importances_
    model_imp_df = pd.DataFrame(pred_imps)
    model_imp_df.index.name = 'Predictor_Index'
    model_imp_df.rename(columns={0: "Relative_Importance"},inplace=True)
    model_imp_df['Predictor_Variable'] = pred_df.columns
    output_dataframe_file = args.prefix + "_Predictor_Importance.tsv"
    model_imp_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    calc_performance_metrics(model=rf, pred_df=pred_df, resp_df=resp_df)
    #MOR = MultiOutputRegressor(rf,n_jobs=1)

def calc_performance_metrics(model=None, pred_df=None, resp_df=None):
    print("Calculating model performance metrics")
    model_info = defaultdict(dict)
    predictions = model.predict(pred_df)
    #print(f"Summary of model predictions:")
    resp_var_names = resp_df.columns
    predictions_df = pd.DataFrame(predictions, columns=resp_var_names)
    predictions_df.index.name = 'Sample_Index'
    predictions_df["Sample"] = pred_df.index
    #predictions_df.index.name = 'Response_Variable'
    output_dataframe_file = args.prefix + "_Model_Training_Set_Predictions.tsv"
    predictions_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    #print(predictions_df.describe())
    for rvar in resp_var_names:
        model_info[rvar]["PearsonR"] = scipy.stats.pearsonr(predictions_df[rvar], resp_df[rvar])[0]
        model_info[rvar]["RMSE"] = np.sqrt(np.mean((predictions_df[rvar] - resp_df[rvar].values.ravel())**2))
    model_info_df = pd.DataFrame.from_dict(model_info)
    model_info_df = model_info_df.transpose()
    model_info_df.index.name = 'Response_Variable'
    output_dataframe_file = args.prefix + "_Model_Performance_Info.tsv"
    model_info_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')
    

@ignore_warnings(category=ConvergenceWarning)
def build_ann(resp_df=None, pred_df=None):
    #Buiild a model for each response variable using the MultiOutputRegressor function
    model_info = defaultdict(dict)
    print(f"Training ANN models")
    #Set up ANN. Hyperparams are the same for all the output variables
    ann = MLPRegressor(hidden_layer_sizes=10, max_iter=1000, learning_rate_init=0.01, random_state=42)
    #regr = ann.fit(pred_df, resp_df)
    #predictions= regr.predict(pred_df)
    # define the direct multioutput MOR model
    MOR = MultiOutputRegressor(ann,n_jobs=args.threads)
    #Calculate the coefficient of determination for each model
    # define the evaluation procedure
    #cv = RepeatedKFold(n_splits=3, n_repeats=1, random_state=1)
    # evaluate the model and collect the scores
    #n_scores = cross_val_score(MOR, pred_df, resp_df, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    # force the scores to be positive
    #n_scores = np.absolute(n_scores)
    # summarize performance
    #print('MAE: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))
    MOR.fit(pred_df, resp_df)
    predictions= MOR.predict(pred_df)
    # make a single prediction
    #row = pred_df.iloc[0]
    #yhat = MOR.predict([row])
    # summarize the prediction
    #print(f"Predicted  1st row {yhat[0]}")
    #print(f"Model performance:")
    #print(f"R^2: {regr.score(sub_res_pre_df[pred_var_names], resp_df)}")
    #print(f"Coefficient of Determination: {regr.score(resp_df, resp_df)}")
    #print(f"Model predictions:")
    #print(predictions)
    print(f"Shape of model predictions: {predictions.shape}")
    resp_var_names = resp_df.columns
    predictions_df = pd.DataFrame(predictions, columns=resp_var_names)
    #print(predictions[1])
    #Calculate the Pearon R2 and RMSE for each response variable
    for rvar in resp_var_names:
        model_info[rvar]["PearsonR"] = scipy.stats.pearsonr(predictions_df[rvar], resp_df[rvar])[0]
        model_info[rvar]["RMSE"] = np.sqrt(np.mean((predictions_df[rvar] - resp_df[rvar].values.ravel())**2))
        # print(f"Model performance for {rvar}:")
        # print(f"\tPearon's R: {scipy.stats.pearsonr(predictions_df[rvar], resp_df[rvar])[0]}")
        # print(f"\tRMSE: {np.sqrt(np.mean((predictions_df[rvar] - resp_df[rvar].values.ravel())**2))}")
    #ann = MLPRegressor(hidden_layer_sizes=10, max_iter=1000, learning_rate_init=0.01, random_state=42)
    #Train ann
    #ann.fit(sub_res_pre_df[pred_var_names], resp_df.values.ravel())
    #Evaluate ann model performance with Peaarsons R squared and RMSE
    #print(f"Model performance:")
    #print(f"R^2: {ann.score(sub_res_pre_df[pred_var_names], resp_df.values.ravel())}")
    model_info_df = pd.DataFrame.from_dict(model_info)
    model_info_df = model_info_df.transpose()
    model_info_df.index.name = 'Response_Variable'
    output_dataframe_file = args.prefix + "_ANN_Niche_Model_Info.tsv"
    model_info_df.to_csv(output_dataframe_file,sep="\t",na_rep='NA')

def format_datasets():
    """Read in the response and predictor tables and format them for model training"""
    print(f"Reading in response data from: {args.response_table}")
    if args.response_index:
        rsp_idx = args.response_index
        print(f"Using: ",rsp_idx," as index column")
    else:
        rsp_idx = 0
        print(f"Using first column as index")

    full_resp_df = pd.read_csv(args.response_table,sep=args.sep_table,index_col=rsp_idx,header=0)
    #print(full_resp_df.describe())

    #Transpose response DF if specified
    if args.transpose_response:
        print("Transposing response table")
        full_resp_df = full_resp_df.transpose()
    #If no response variables are specified, use all variables in the table
    resp_var_names = []
    if args.response_vars is None:
        print("Using all response variables")
        resp_var_names = full_resp_df.columns
    else:    
        print(f"Using response variables: {args.response_vars}")
        resp_var_names = args.response_vars

    sub_resp_df = full_resp_df[resp_var_names]

    print(f"Reading in predictor data from: {args.predictor_table}")
    if args.predictor_index:
        prd_idx = args.predictor_index
        print(f"Using: ",prd_idx," as index column")
    else:
        prd_idx = 0
        print(f"Using first column as index")

    full_pred_df = pd.read_csv(args.predictor_table,sep=args.sep_table,index_col=prd_idx,header=0)
    #print(full_pred_df.describe())

    #Transpose predictor DF if specified
    if args.transpose_predictors:
        print("Transposing predictor table")
        full_pred_df = full_pred_df.transpose()

    #If no predictor variables are specified, use all variables in the table
    if args.predictor_vars is None:
        print("Using all predictor variables")
        pred_var_names = full_pred_df.columns
    else:
        print(f"Using predictor variables: {args.predictor_vars}")
        pred_var_names = args.predictor_vars

    sub_pred_df = full_pred_df[pred_var_names]

    #Concatenate predictor and response variables
    print("Concatenating predictor and response variables")
    sub_res_pre_df = pd.concat([sub_resp_df, sub_pred_df], axis=1)
    print(f"Shape of DF after contatenation: {sub_res_pre_df.shape}")
    #print(sub_res_pre_df.describe())
    #Remove rows with missing values in sub_res_pre_df
    print("Removing rows with missing values")
    sub_res_pre_df.dropna(inplace=True)
    #Print shape of sub_res_pre_df
    print(f"Shape of DF passed to model training: {sub_res_pre_df.shape}")
    resp_df = sub_res_pre_df[resp_var_names]
    pred_df = sub_res_pre_df[pred_var_names]
    return(resp_df, pred_df)

central()
