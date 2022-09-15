from __future__ import print_function
import IMP
import GoodScoringModelSelector
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="List and extract good-scoring models from a set of sampling runs. Example of usage: select_good_scoring_models.py -rd <run_directory_for_sampling> -rp <run_prefix> -sl ExcludedVolumeSphere_None GaussianEMRestraint_None -pl CrossLinkingMassSpectrometryDataScore|XLDSS CrossLinkingMassSpectrometryDataScore|XLEDC -agl -9999999.0 -99999.0 -aul 99999999.0 999999.0 -mlt 0 0 -mut 0 0. Flag -h for more details.")
    
    parser.add_argument("-rd","--run_directory",dest="run_dir",help="directory in which sampling results are stored") 
    
    parser.add_argument("-rp","--run_prefix",dest="run_prefix",help="prefix of runs") 
                        
    parser.add_argument("-sl","--selection_keywords_list",nargs='+',type=str,dest="selection_keywords_list",help="list of stat file keywords corresponding to selection criteria")
    parser.add_argument("-pl","--printing_keywords_list",nargs='+',type=str,dest="printing_keywords_list",help="list of stat file keywords whose values are printed out for selected models")
    
    # thresholds only apply to selection keywords
    parser.add_argument("-alt","--aggregate_lower_thresholds",nargs='+',type=float,dest="aggregate_lower_thresholds",help="aggregate lower thresholds")
    parser.add_argument("-aut","--aggregate_upper_thresholds",nargs='+',type=float,dest="aggregate_upper_thresholds",help="aggregate upper thresholds")
    parser.add_argument("-mlt","--member_lower_thresholds",nargs='+',type=float,dest="member_lower_thresholds",help="member lower thresholds")
    parser.add_argument("-mut","--member_upper_thresholds",nargs='+',type=float,dest="member_upper_thresholds",help="member upper thresholds")

    parser.add_argument("-e","--extract",default=False,dest="extract",action='store_true',help="Type -e to extract all good scoring model RMFs from the trajectory files")
    parser.add_argument("-sf","--score_file",default="scores", type=str, dest="score_file_prefix",help="Score file prefix for samples A and B. Default is 'scores'")
    result = parser.parse_args()
 
    return result
    
def select_good_scoring_models():
     
    # process input
    arg=parse_args()
    print(arg.selection_keywords_list)
    gsms=GoodScoringModelSelector.GoodScoringModelSelector(arg.run_dir,arg.run_prefix)
               
    subsets = gsms.get_good_scoring_models(selection_keywords_list=arg.selection_keywords_list,printing_keywords_list=arg.printing_keywords_list,
    aggregate_lower_thresholds=arg.aggregate_lower_thresholds,aggregate_upper_thresholds=arg.aggregate_upper_thresholds,
    member_lower_thresholds=arg.member_lower_thresholds,member_upper_thresholds=arg.member_upper_thresholds,extract=arg.extract)
    return subsets
        
def create_score_files(subsets, field="Total_Score"):
    arg=parse_args()
    scoreA = open("good_scoring_models/" + arg.score_file_prefix + "A.txt","w")
    scoreB = open("good_scoring_models/" + arg.score_file_prefix + "B.txt","w")

    if not arg.extract:
        model_file = open("filter/model_ids_scores.txt","r")
    else:
        model_file = open("good_scoring_models/model_ids_scores.txt","r")
    print("Creating input files for Total_Score convergence test")

    for line_index,each_model_line in enumerate(model_file.readlines()):

        # Find index of the field we want to use for model score convergence
        if line_index==0:
            field_headers = each_model_line.strip().split()
            ts_ix = field_headers.index(field)
            run_ix = field_headers.index("Run_id")
            model_ix = field_headers.index("Model_index")


        else:
            fields = each_model_line.strip().split()
            score=fields[ts_ix]
            if arg.extract:
                model = int(fields[model_ix])
                print(score, file=scoreA if model in subsets[0] else scoreB)
            else:
                if int(fields[run_ix])==1:
                    print(score, file=scoreA)
                elif int(fields[run_ix])==2:
                    print(score, file=scoreB)
                else:
                    print("create_scores_file: model_ids_scores.txt file has an incorrect format.")
                    exit()
    scoreA.close()
    scoreB.close()
    

    
 
if __name__ == "__main__" :
    
    arg=parse_args()
    subsets = select_good_scoring_models()

    # Create Score Files
    if arg.extract:
        create_score_files(subsets)

        print("Ready to calculate sampling precision with Master_Sampling_Exhaustiveness_Analysis.py")
    else:
        print("Model score file at ./filter/model_ids_scores.txt")
