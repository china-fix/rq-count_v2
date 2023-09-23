#!/usr/bin/env python3


'''
This script is design for caculate the relative mapping rate change between start samples and end samples
### 20220512-update modify the linear model, former (block depth)X > (block depth)Y, updated (block depth, block position)X > (block depth)Y
### 20220516-update add new machine-learning model to instead the former linear model
### 20220601-update major revision, the prediction strategy has modified!
### 20220604-update update the outlier dropping function!
### 20230507-update add multiprocessing to use mutiple CPU cores; handle BLAST no match issues; handle other unexpected errors
'''
# import glob
# import os
import subprocess
# import copy
import sys
import argparse
import io
from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
# from sklearn import linear_model
import numpy as np

import multiprocessing
from functools import partial

# from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
# from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import r2_score, explained_variance_score, mean_absolute_error, mean_absolute_percentage_error, median_absolute_error
# from sklearn.preprocessing import PolynomialFeatures
from sklearn.svm import SVR


# from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--MAPPING_REF', required=True, type=str, metavar='FILENAME', help="the fasta filename you used for bwa mapping reference")
    parser.add_argument('--TARGETS', required=True, type=str, metavar='FILENAME', help="the fasta file includes all the tested gene sequences you want to measure")
    parser.add_argument('--DEPTH_TAB', required=True, type=str, metavar='FILENAME', help="the tab file describe the mapping depth of every base")
    # parser.add_argument('--DEPTH_TAB_E', required=True, type=str, metavar='FILENAME', help="the tab file (end condiction) describe the mapping depth of every base")
    parser.add_argument('--CUTLEN', default=5000, type=int, metavar='DEFAULT 5000', help="the flanking seq length for target catulation part (extracted for model learning)")
    parser.add_argument('--FIXLEN', type=int, help="flanking seq length you want to drop near the deletion edge location (default is 0)")
    # parser.add_argument('--DROP_ABOVE', default=1000, type=int, metavar='DEFAULT 1000', help="mapping depth above this number will be dropped")
    # parser.add_argument('--DROP_BELOW', default=0, type=int, metavar='DEFAULT 0', help="mapping depth below this number will be dropped")

    parser.add_argument('--OUT', default="report", type=str, metavar='FILENAME', help="Output file name, report as default")
    parser.add_argument('--DEV', action='store_const', const=True, metavar='FOR DEVELOPPING ONLY', help="for developping and testing only normal user can just ignore this")
    # parser.add_argument('--TEST_MODE', action='store_const', const=True, metavar='FOR TEST', help="only for coding test")
    # parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    return parser.parse_args()



### get targets locations, mask these locations and output df for scikit-learn fit
### same time also prepare the masked locations and output df for scikit-learn predict later
def get_learning_df(depth_tab, targets, mapping_ref, cutlen, fixlen=False):
    
    depth_tab_df = pd.read_csv(depth_tab, sep='\t', names=["reference","1_index","depth"])

    ## doing blast to get the xml_obj
    blast_run_xml = subprocess.run(["blastn", "-query", targets, "-subject", mapping_ref, "-outfmt", "5"], check=True, capture_output=True)
    xml = blast_run_xml.stdout
    xml_obj = io.BytesIO(xml)
    # xml_obj_c = io.BytesIO(xml)

    ## read xml_obj and get the locations and remove these parts in depth_tab_df, same time get the masked_df
    blast_records = NCBIXML.parse(xml_obj)

    masked_df_append = pd.DataFrame()
    flanking_df_append = pd.DataFrame()

    for blast_record in blast_records:
        try:
            alignment = blast_record.alignments[0]
        except IndexError:
            print(blast_record.query +"BLAST no match!")
            continue
        # print(blast_record.query +"BLAST match!")
        hsp = alignment.hsps[0]
        start_loc = hsp.sbjct_start   ## here followed 1-index
        end_loc = hsp.sbjct_end
        target_name = blast_record.query
        if start_loc > end_loc:
            start_loc, end_loc = end_loc, start_loc
        if fixlen:
            print("Hi Xiao i will use the fixlen mode to prepare the data, the fixlen is: " + str(fixlen))
            masked_df = depth_tab_df.query('(@start_loc) <= `1_index` <= (@end_loc)')
            masked_df = masked_df.assign(target_name=target_name)
            masked_df_append = pd.concat([masked_df_append, masked_df])
            # depth_tab_df = depth_tab_df.query('not (@start_loc + @fixlen) <= `1_index` <= (@end_loc-@fixlen)')
            # flanking_df = depth_tab_df.query('((@end_loc + @cutlen) >= `1_index` > @end_loc) or ((@start_loc - @cutlen) <= `1_index` < @start_loc)')
            flanking_df = depth_tab_df.query('((@end_loc + @cutlen) >= `1_index` > (@end_loc+@fixlen)) or ((@start_loc - @cutlen) <= `1_index` < (@start_loc - @fixlen))')
            flanking_df = flanking_df.assign(target_name=target_name)
            flanking_df_append = pd.concat([flanking_df_append, flanking_df]) 
        else:
            print("Hi Xiao will not use fixlen mode to prepare the data")
            masked_df = depth_tab_df.query('@start_loc <= `1_index` <= @end_loc')
            masked_df = masked_df.assign(target_name=target_name)
            masked_df_append = pd.concat([masked_df_append, masked_df])
            # depth_tab_df = depth_tab_df.query('not @start_loc <= `1_index` <= @end_loc')
            flanking_df = depth_tab_df.query('((@end_loc + @cutlen) >= `1_index` > @end_loc) or ((@start_loc - @cutlen) <= `1_index` < @start_loc)')
            flanking_df = flanking_df.assign(target_name=target_name)
            flanking_df_append = pd.concat([flanking_df_append, flanking_df]) 
    
    
    # ## read xml_obj_c and get the locations and get the flanking_df
    # blast_records = NCBIXML.parse(xml_obj_c)

    # flanking_df_append = pd.DataFrame()
    # for blast_record in blast_records:
    #     try:
    #         alignment = blast_record.alignments[0]
    #     except IndexError:
    #         continue
    #     hsp = alignment.hsps[0]
    #     start_loc = hsp.sbjct_start   ## here followed 1-index
    #     end_loc = hsp.sbjct_end
    #     target_name = blast_record.query
    #     if start_loc > end_loc:
    #         start_loc, end_loc = end_loc, start_loc

    #     flanking_df = depth_tab_df.query('((@end_loc + @cutlen) >= `1_index` > @end_loc) or ((@start_loc - @cutlen) <= `1_index` < @start_loc)')
    #     # flanking_df['target_name'] = target_name
    #     flanking_df = flanking_df.assign(target_name=target_name)
    #     flanking_df_append = pd.concat([flanking_df_append, flanking_df])      

    return depth_tab_df, masked_df_append, flanking_df_append


### get the ML model
def get_ML_model (learning_df):
    X = learning_df[["1_index"]].to_numpy(dtype=float)
    y = learning_df["depth"].to_numpy(dtype=float)


    regr = SVR(kernel="rbf", C=50, gamma= 0.1, epsilon=0.1, cache_size=5000)
    # regr = SVR(kernel="linear", C=100, gamma="auto", epsilon=0.1, cache_size=4000)
    # regr = MLPRegressor(random_state=1, max_iter=2000, early_stopping=True, solver="adam")
    regr.fit(X, y)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)
    regr_metrics = SVR(kernel="rbf", C=50, gamma= 0.1, epsilon=0.1, cache_size=5000)
    # regr_metrics = SVR(kernel="linear", C=100, gamma="auto", epsilon=0.1, cache_size=4000)
    # regr_metrics = MLPRegressor(random_state=1, max_iter=2000, early_stopping=True,solver="adam")
    
    regr_metrics.fit(X_train, y_train)
    R2 = r2_score(y_test, regr_metrics.predict(X_test))
    print("R2:", R2)
    VAR = explained_variance_score(y_test, regr_metrics.predict(X_test))
    print("variance:", VAR)
    ABS = mean_absolute_error(y_test, regr_metrics.predict(X_test))
    print("mean absolute error:", ABS)
    ABS_MEDIAN = median_absolute_error(y_test, regr_metrics.predict(X_test))
    print("median absolute error:", ABS_MEDIAN)
    ABS_PER = mean_absolute_percentage_error(y_test, regr_metrics.predict(X_test))
    print("mean absolute percentage error: ", ABS_PER)

    return regr, ABS, R2, ABS_MEDIAN

### input pandas series and output series without outliters
def elimate_outliters (data):
    mu = data.mean()
    std = data.std()
    # sem = std/np.sqrt(len(data))
    clean_data = data [data > mu - 1.645 * std] [data < mu + 1.645*std]
    # clean_data = data [data > mu - 1.645 * sem] [data < mu + 1.645*sem]

    return clean_data



### use trained model to predict, then caculate the ratio
def ratio_caculate (model, targets_df, mean_abs_er, R2, median_abs_er):
    
    predict_X = targets_df[[ "1_index"]].to_numpy(dtype=float)

    # compare_Y = df_after_targets["depth_sum"].to_numpy(dtype=float)
    # predict_X_reshape = predict_X.reshape(-1, 1)
    predict_y = model.predict(predict_X)
    # targets_df_new = pd.concat([targets_df, pd.DataFrame(predict_y, columns=["predict_depth"])])    
    # targets_df['predict_depth'] = predict_y.tolist()
    targets_df = targets_df.assign(predict_depth=predict_y.tolist())
    targets_df['DEL_depth'] = targets_df.predict_depth - targets_df.depth
    targets_df['target_rate'] = (targets_df.predict_depth-targets_df.depth)/targets_df.predict_depth
    

    #doing the summary
    clean_df = targets_df[["target_name", "depth", "predict_depth","DEL_depth","target_rate"]]


    # df=pd.DataFrame({"target_name":[0], "depth":[0], "predict_depth":[0]})
    df=clean_df.groupby("target_name").median()
    
    ## re-caculate the median depth and median predict_depth
    series_depth = clean_df["depth"]
    series_predict_depth = clean_df["predict_depth"]
    series_DEL_depth = clean_df["DEL_depth"]
    series_target_rate =clean_df["target_rate"]
    # cleaning
    series_depth_clean = elimate_outliters(series_depth)
    series_predict_depth_clean = elimate_outliters(series_predict_depth)
    series_DEL_depth_clean = elimate_outliters(series_DEL_depth)
    series_target_rate_clean = elimate_outliters(series_target_rate)
    # assigning
    df = df.assign(depth_c=series_depth_clean.median())
    df = df.assign(depth_c_std=series_depth_clean.std())
    df = df.assign(predict_depth_c=series_predict_depth_clean.median())
    df = df.assign(predict_depth_c_std=series_predict_depth_clean.std())
    df = df.assign(DEL_depth_c=series_DEL_depth_clean.median())
    df = df.assign(DEL_depth_c_std=series_DEL_depth_clean.std())
    # df = df.assign(target_rate_c_sem=series_target_rate_clean.std() / np.sqrt(len(series_target_rate_clean)))
    df = df.assign(target_rate_c=series_target_rate_clean.median())
    df = df.assign(target_rate_c_std=series_target_rate_clean.std())

    # df = df.assign(target_rate_c_sem=series_target_rate_clean.std / np.sqrt(len(series_target_rate_clean)))
    # df['DEL_depth'] = df.predict_depth - df.depth
    ##assign Machine learning evaluation number
    df = df.assign(ML_median_abs_er=median_abs_er)
    df = df.assign(ML_mean_abs_er=mean_abs_er)
    df = df.assign(ML_R2=R2)
    # df['target_rate'] = (df.predict_depth-df.depth)/df.predict_depth
    


    return targets_df, df

 

def process_target(name, flanking_df, targets_df):
    try:
        print("\nHi xiao, caculation for " + name)
        # get the model
        flanking_df_i = flanking_df.query('`target_name` == @name')
        trained_model, mean_abs_er, R2, median_abs_er = get_ML_model(flanking_df_i)

        ## use linear/ML model to predict targets depth, then caculate
        targets_df_i = targets_df.query('`target_name` == @name')
        new_df, new_df_summary =ratio_caculate(trained_model, targets_df_i, mean_abs_er, R2, median_abs_er)

        ## output the flanking_df_i plus target_df_i
        flanking_df_i['target_name']='flanking_'+name
        new_df=pd.concat([new_df,flanking_df_i])


        return new_df, new_df_summary
    
    except Exception as e:
        print(f"Error processing target {name}: {e}")
        return None, None


def main():
    args = parse_args()
    new_df_append = pd.DataFrame()
    new_df_append_s = pd.DataFrame()
    if args.DEV:
        print("Hi Xiao, i am in testing mode!!!")
        #prepare the feeding data
        feeding_df, targets_df, flanking_df = get_learning_df(args.DEPTH_TAB, args.TARGETS, args.MAPPING_REF, args.CUTLEN, args.FIXLEN)
        # get the model
        trained_model = get_ML_model(feeding_df)
        ## use linear/ML model to predict targets depth, then caculate
        new_df_append=ratio_caculate(trained_model, targets_df)
    else:
        #prepare the feeding data
        feeding_df, targets_df, flanking_df = get_learning_df(args.DEPTH_TAB, args.TARGETS, args.MAPPING_REF, args.CUTLEN, args.FIXLEN)
        print("Hi Xiao, data preparation complete. ")

        #learning the flanking parts followed with prediction individually for every targets 
        # for name in flanking_df.target_name.unique():
        #     print("\nHi xiao, caculation for " + name)
        #     # get the model
        #     flanking_df_i = flanking_df.query('`target_name` == @name')
        #     trained_model, mean_abs_er, R2, median_abs_er = get_ML_model(flanking_df_i)

        #     ## use linear/ML model to predict targets depth, then caculate
        #     targets_df_i = targets_df.query('`target_name` == @name')
        #     new_df, new_df_summary =ratio_caculate(trained_model, targets_df_i, mean_abs_er, R2, median_abs_er)
        #     new_df_append = pd.concat([new_df_append, new_df])
        #     new_df_append_s = pd.concat([new_df_append_s, new_df_summary])

       

       
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()-1) as pool:
        partial_func = partial(process_target, flanking_df=flanking_df, targets_df=targets_df)
        results = pool.map(partial_func, flanking_df.target_name.unique())

        new_df_append = pd.concat([result[0] for result in results])
        new_df_append_s = pd.concat([result[1] for result in results])
                    
            
    print("\nHi Xiao, I am done and i will save the report!")
    new_df_append.to_csv(args.OUT+"_detailed.csv")
    new_df_append_s.to_csv(args.OUT+"_summary.csv")
    print("detailed report file named: " + args.OUT+"_detailed.csv")
    print("summary report file named: " + args.OUT+"_summary.csv")
    
    





if __name__ == '__main__':
    sys.exit(main())