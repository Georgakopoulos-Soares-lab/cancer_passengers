import os
import argparse
import pandas as pd 

# input files/dir
FABIAN_OUTPUT = "data/fabian_output"

# output files/dir
TFBS_DATA_DIFF = "data/TFBS_data_diff"
TFBS_DATA_SCORES = "data/TFBS_data_scores"

if not os.path.exists(TFBS_DATA_DIFF):
    os.makedirs(TFBS_DATA_DIFF)

if not os.path.exists(TFBS_DATA_SCORES):
    os.makedirs(TFBS_DATA_SCORES)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer_type", type=str, help="Cancer tissue type e.g Pancreas")
    args = parser.parse_args()
    cancer_type = args.cancer_type
    print(f"Processing FABIAN output for {cancer_type}...")
    
    # FABIAN output files for the cancer type
    files = os.listdir(FABIAN_OUTPUT)
    files = [f for f in files if cancer_type in f]
    print(files)

    # process FABIAN output
    fab_out = pd.DataFrame()
    for f in files:
        df = pd.read_csv(f"{FABIAN_OUTPUT}/{f}", sep="\t")
        headers = ["variant", "tf", "model_id", "model_db", "wt_score", "mt_score", "start_wt", "end_wt", "start_mt", "end_mt", "strand_wt", "strand_mt", "prediction", "score"]
        df.columns = headers
        print(df.shape, fab_out.shape)
        fab_out = pd.concat([fab_out, df], axis=0)
    fab_out["variant"] = fab_out["variant"].apply(lambda x: x.split(".")[0])
    fab_out.drop_duplicates(inplace=True)

    # mean difference in TF binding effect of WT and mutant sequences
    fab_out_diff = fab_out.groupby(["variant", "tf"]).agg({
        "score": "mean"
    }).reset_index()
    fab_out_diff = fab_out_diff[fab_out_diff["score"] != 0]
    fab_out_diff.to_csv(f"{TFBS_DATA_DIFF}/{cancer_type}.tsv", sep="\t", index=False)

    # mean WT and mutant scores 
    fab_out_scores = fab_out.groupby(["variant", "tf"]).agg({
        "wt_score": "mean",
        "mt_score": "mean",
    }).reset_index()
    fab_out_scores["diff_score"] = fab_out_scores["mt_score"] - fab_out_scores["wt_score"]
    fab_out_scores.to_csv(f"{TFBS_DATA_SCORES}/{cancer_type}.tsv", sep="\t", index=False)
