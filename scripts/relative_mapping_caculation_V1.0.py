#!/usr/bin/env python3


'''
This script is design for caculate the relative mapping rate change between start samples and end samples
### 20220512-update modify the linear model, former (block depth)X > (block depth)Y, updated (block depth, block position)X > (block depth)Y
### 20220516-update add new machine-learning model to instead the former linear model
'''
##old version no longer support
# import glob
import os
import subprocess
# import copy
import sys
import argparse
import io
from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
from sklearn import linear_model
import numpy as np



from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import r2_score, explained_variance_score, mean_absolute_error, mean_absolute_percentage_error
from sklearn.preprocessing import PolynomialFeatures


# from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--MAPPING_REF', required=True, type=str, metavar='FILENAME', help="the fasta filename you used for bwa mapping reference")
    parser.add_argument('--TARGETS', required=True, type=str, metavar='FILENAME', help="the fasta file includes all the tested gene sequences you want to measure")
    parser.add_argument('--DEPTH_TAB_S', required=True, type=str, metavar='FILENAME', help="the tab file (start condiction) describe the mapping depth of every base")
    parser.add_argument('--DEPTH_TAB_E', required=True, type=str, metavar='FILENAME', help="the tab file (end condiction) describe the mapping depth of every base")
    parser.add_argument('--CUTLEN', default=200, type=int, metavar='DEFAULT 200', help="the separation interval length for the MAPPING_REF")

    parser.add_argument('--DROP_ABOVE', default=1000, type=int, metavar='DEFAULT 1000', help="mapping depth above this number will be dropped")
    parser.add_argument('--DROP_BELOW', default=0, type=int, metavar='DEFAULT 0', help="mapping depth below this number will be dropped")

    parser.add_argument('--OUT', default="fold_change.csv", type=str, metavar='FILENAME', help="Output file name, fould_change.csv as default")
    parser.add_argument('--ML', action='store_const', const=True, metavar='MACHINE LEARNING', help="use machine learning model")
    parser.add_argument('--TEST_MODE', action='store_const', const=True, metavar='FOR TEST', help="only for coding test")
    # parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    return parser.parse_args()


### make ref.bed
def make_ref_bed(mapping_ref, cutlen):
    ## get ref length
    mapping_ref_obj = SeqIO.read(mapping_ref, "fasta")
    range_end = len(mapping_ref_obj.seq)
    start_loc = 0
    for end_loc in range(cutlen, range_end, cutlen):
        print("ref\t"+str(end_loc-cutlen)+"\t"+str(end_loc)+"\tBLOCK_"+str(end_loc)+"\t"+str((start_loc +end_loc)/2), file=open("temp_ref.bed", "a") )
        start_loc=end_loc


### make targets.bed
def make_targets_bed(targets, mapping_ref):
    ## doing blast to get the xml_obj
    blast_run_xml = subprocess.run(["blastn", "-query", targets, "-subject", mapping_ref, "-outfmt", "5"], check=True, capture_output=True)
    xml = blast_run_xml.stdout
    xml_obj = io.BytesIO(xml)
    
    ## read xml_obj and get the locations
    blast_records = NCBIXML.parse(xml_obj)
    for blast_record in blast_records:
        alignment = blast_record.alignments[0]
        hsp = alignment.hsps[0]
        start_loc = hsp.sbjct_start-1
        end_loc = hsp.sbjct_end
        target_name = blast_record.query
        if start_loc < end_loc:
            print("ref\t"+str(start_loc)+"\t"+str(end_loc)+"\tTARGET_"+target_name+"\t"+str((start_loc +end_loc)/2), file=open("temp_targets.bed", "a") )
        else:
            print("ref\t"+str(end_loc)+"\t"+str(start_loc)+"\tTARGET_"+target_name+"\t"+str((start_loc +end_loc)/2), file=open("temp_targets.bed", "a") )


### read bed locations and output df for scikit-learn fit
def sum_block_depth(bedfile, depth_tab):
    depth_tab_df = pd.read_csv(depth_tab, sep='\t', names=["reference","1_index","depth"])
    bedfile_df = pd.read_csv(bedfile, sep='\t', names=["reference","start","end","block","block_index"])
    
    rows = []
    for ind in bedfile_df.index:
        start_num = bedfile_df['start'][ind]
        end_num = bedfile_df['end'][ind]
        block_depth_sum = depth_tab_df.loc[start_num : end_num-1]['depth'].mean()   ###change to mean

        block_name = bedfile_df['block'][ind]
        block_index = bedfile_df['block_index'][ind]

        rows.append([block_name, block_depth_sum, block_index])
    
    fit_df = pd.DataFrame(rows, columns=["block_name", "depth_sum", "block_index"])
    # print(fit_df)
    return fit_df

### get the linear model
def get_linear_model (df_before, df_after):
    linear_model_X  = df_before[["depth_sum","block_index"]].to_numpy(dtype=float)
    # linear_model_X = linear_model_X.reshape(-1, 1)
    # print(linear_model_X)
    linear_model_Y  = df_after["depth_sum"].to_numpy(dtype=float)
    # print(linear_model_Y)
    
    # Create linear regression object
    regr = linear_model.LinearRegression()
    # Train the model using sets
    regr.fit(linear_model_X, linear_model_Y)

    ## model metrics
    X_train, X_test, y_train, y_test = train_test_split(linear_model_X, linear_model_Y, random_state=1)
    regr_metrics = linear_model.LinearRegression()
    regr_metrics.fit(X_train, y_train)
    R2 = r2_score(y_test, regr.predict(X_test))
    print(R2)
    
    return regr

### get the ML model
def get_ML_model (df_before, df_after):
    ML_model_X  = df_before[["depth_sum","block_index"]].to_numpy(dtype=float)
    
    ML_model_Y  = df_after["depth_sum"].to_numpy(dtype=float)
    # print(linear_model_Y)
    
    # Create ML regression object
    # regr = MLPRegressor(random_state=1, max_iter=500)
    regr = KernelRidge()

    # Train the model using sets
    regr.fit(ML_model_X, ML_model_Y)

     ## model metrics
    X_train, X_test, y_train, y_test = train_test_split(ML_model_X, ML_model_Y, random_state=1)
    regr_metrics = KernelRidge()
    regr_metrics.fit(X_train, y_train)
    R2 = r2_score(y_test, regr.predict(X_test))
    print(R2)

    return regr

### get the TEST model
def get_TEST_model (df_before, df_after, drop_below=0, drop_above=1000):
    # add pre-data-processing to make the learning 
    df_merge = df_before.merge(df_after, on='block_index')
    # print(df_merge)
    df_merge = df_merge.query('( @drop_above> depth_sum_x > @drop_below) and (@drop_above > depth_sum_y > @drop_below)')
    # print(df_merge)
    # ML_model_X  = df_merge[["depth_sum_x"]].to_numpy(dtype=float)
    ML_model_X  = df_merge[["depth_sum_x","block_index"]].to_numpy(dtype=float)
    # print(ML_model_X)
    
    ML_model_Y  = df_merge["depth_sum_y"].to_numpy(dtype=float)
    # print(ML_model_Y)

    # ML_model_X  = df_before[["depth_sum","block_index"]].to_numpy(dtype=float)
    
    # ML_model_Y  = df_after["depth_sum"].to_numpy(dtype=float)
    # print(linear_model_Y)
    
    # Create ML regression object
    # regr = linear_model.LinearRegression()
    # regr = MLPRegressor(random_state=1, max_iter=2000, early_stopping=True, hidden_layer_sizes=[100,100],solver="lbfgs")
    from sklearn.svm import SVR
    regr = SVR(kernel="rbf", C=1, gamma="auto", degree=3, epsilon=0.1, coef0=0, cache_size=3000)
    
    # regr = KernelRidge()

    # Train the model using sets
    regr.fit(ML_model_X, ML_model_Y)

     ## model metrics
    X_train, X_test, y_train, y_test = train_test_split(ML_model_X, ML_model_Y, random_state=1)
    # print(y_test)
    # print(regr.predict(X_test))
    # regr_metrics = MLPRegressor(random_state=1, max_iter=2000, early_stopping=True, hidden_layer_sizes=[100,100],solver="lbfgs")
    regr_metrics = SVR(kernel="rbf", C=1, gamma="auto", degree=3, epsilon=0.1, coef0=0, cache_size=3000)
    
    regr_metrics.fit(X_train, y_train)
    R2 = r2_score(y_test, regr_metrics.predict(X_test))
    print(R2)
    VAR = explained_variance_score(y_test, regr_metrics.predict(X_test))
    print(VAR)
    ABS = mean_absolute_error(y_test, regr_metrics.predict(X_test))
    print(ABS)
    ABS_PER = mean_absolute_percentage_error(y_test, regr_metrics.predict(X_test))
    print(ABS_PER)
    
    # temp=y_test/regr_metrics.predict(X_test)
    # for row in temp:
    #     print(row, file=open("temp_fold", "a") )
   

    
    # temp=ML_model_Y/regr.predict(ML_model_X)
    # for row in temp:
    #     print(row, file=open("temp_fold", "a") )
   

    return regr

### get model
def get_model (df_before, df_after, drop_below=0, drop_above=1000, ml=False):
    # add pre-data-processing to make the learning 
    df_merge = df_before.merge(df_after, on='block_index')
    # print(df_merge)
    df_merge = df_merge.query('( @drop_above> depth_sum_x > @drop_below) and (@drop_above > depth_sum_y > @drop_below)')
    # print(df_merge)
    # ML_model_X  = df_merge[["depth_sum_x"]].to_numpy(dtype=float)
    model_X  = df_merge[["depth_sum_x","block_index"]].to_numpy(dtype=float)
    # print(ML_model_X)
    
    model_Y  = df_merge["depth_sum_y"].to_numpy(dtype=float)
    # print(ML_model_Y)
    
    # Create ML regression object
    if ml:
        regr = MLPRegressor(random_state=1, max_iter=2000, early_stopping=True,solver="adam")
        print("hi Xiao i will use machine learning model to learn")
        ## model metrics
        X_train, X_test, y_train, y_test = train_test_split(model_X, model_Y, random_state=1)
        regr_metrics = MLPRegressor(random_state=1, max_iter=20000, early_stopping=True,solver="adam")
        regr_metrics.fit(X_train, y_train)
        R2 = r2_score(y_test, regr_metrics.predict(X_test))
        print("R2 score is:", R2)
        VAR = explained_variance_score(y_test, regr_metrics.predict(X_test))
        print("variance score is:", VAR)
        ABS = mean_absolute_error(y_test, regr_metrics.predict(X_test))
        print("mean absolute error score is:", ABS)
        ABS_PER = mean_absolute_percentage_error(y_test, regr_metrics.predict(X_test))
        print("mean absolute percentage error is: ", ABS_PER)
    else:
        regr = linear_model.LinearRegression()
        print("hi Xiao i will use linear model to learn")
        ## model metrics
        X_train, X_test, y_train, y_test = train_test_split(model_X, model_Y, random_state=1)
        regr_metrics = linear_model.LinearRegression()
        regr_metrics.fit(X_train, y_train)
        R2 = r2_score(y_test, regr_metrics.predict(X_test))
        print(R2)
        VAR = explained_variance_score(y_test, regr_metrics.predict(X_test))
        print(VAR)
        ABS = mean_absolute_error(y_test, regr_metrics.predict(X_test))
        print(ABS)
        ABS_PER = mean_absolute_percentage_error(y_test, regr_metrics.predict(X_test))
        print(ABS_PER)
    
    # Train the model using sets
    regr.fit(model_X, model_Y)
    return regr
    # temp=ML_model_Y/regr.predict(ML_model_X)
    # for row in temp:
    #     print(row, file=open("temp_fold", "a") )
   

    # return regr



### use linear/ML model to predict after level, then caculate the fold change
def fold_change_caculate (model, df_before_targets, df_after_targets):
    targets_name_s = df_before_targets["block_name"] #Series
    
    predict_X = df_before_targets[["depth_sum", "block_index"]].to_numpy(dtype=float)
    compare_Y = df_after_targets["depth_sum"].to_numpy(dtype=float)
    # predict_X_reshape = predict_X.reshape(-1, 1)
    predict_Y = model.predict(predict_X)
    fold_change_Y = compare_Y/predict_Y


    # predict_X_s = pd.Series(predict_X, name="before_depth_sum")
    compare_Y_s = pd.Series(compare_Y, name="after_depth")
    predict_Y_s = pd.Series(predict_Y, name="predicted_after_depth")
    fold_change_Y_s = pd.Series(fold_change_Y, name="fold_change")

    df = pd.concat([targets_name_s,predict_Y_s,compare_Y_s,fold_change_Y_s], axis=1)
    return df
    




def main():
    args = parse_args()
    ## prepare ref.bed
    make_ref_bed(args.MAPPING_REF, args.CUTLEN)
    print("Hi Xiao, i made a temp_ref.bed file")

    ## prepare targets.bed
    make_targets_bed(args.TARGETS, args.MAPPING_REF)
    print("Hi Xiao, i made a temp_targets.bed file")

    ## runing bedtools, bedtools intersect -a ref.bed -b targets.bed -v
    bedtools_run = subprocess.run(["bedtools", "intersect", "-a", "temp_ref.bed", "-b", "temp_targets.bed", "-v"], check=True, capture_output=True)
    bed_obj = io.BytesIO(bedtools_run.stdout)
    bed_obj_c = io.BytesIO(bedtools_run.stdout)
    # print(bed_obj)
    ## get block depth sum of before and after, then feed to liner model
    feed_df_after = sum_block_depth(bed_obj, args.DEPTH_TAB_E)
    feed_df_before = sum_block_depth(bed_obj_c, args.DEPTH_TAB_S)
    

    # get the model or testing
    if args.TEST_MODE:
        trained_model = get_TEST_model(feed_df_before, feed_df_after, args.DROP_BELOW, args.DROP_ABOVE)
        print("Hi Xiao, I learned TEST model from your feed data")
    else:
        trained_model = get_model(feed_df_before, feed_df_after, args.DROP_BELOW, args.DROP_ABOVE, args.ML)
        print("Hi Xiao, I learned model from your feed data")


    ## use linear/ML model to predict targets after level, then caculate the fold change
    caculate_df_before = sum_block_depth("temp_targets.bed", args.DEPTH_TAB_S)
    caculate_df_after = sum_block_depth("temp_targets.bed", args.DEPTH_TAB_E)
    
    
    fold_change_df = fold_change_caculate(trained_model, caculate_df_before, caculate_df_after)
    print("Hi Xiao, I am done and i will save the caculation file named as "+ args.OUT)
    fold_change_df.to_csv(args.OUT)
    os.system("mkdir --parents ./temp_bed; mv temp_targets.bed ./temp_bed/. ; mv temp_ref.bed ./temp_bed/.")
    print("Ohh, I also move the temp files to temp_bed folder")
    





if __name__ == '__main__':
    sys.exit(main())