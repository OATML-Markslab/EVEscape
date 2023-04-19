import pandas as pd
import numpy as np
import math
import argparse

# Initiate the parser
parser = argparse.ArgumentParser()

# Add long and short argument
parser.add_argument("--evescape", "-e", help="Path to EVEscape scores")
parser.add_argument("--gisaid", "-g", help="Path to GISAID summary table")

# Read arguments from the command line
args = parser.parse_args()

EVEscape_path = args.evescape
gisaid_path = args.gisaid


def sigmoid(x):
    return 1 / (1 + math.exp(-x))


def clean_up_mutations(mut_list):

    indel_mapper = {
        "LPPA24S": "A27S",
        "EFR156G": "E156G",
        "GVYY142D": "G142D",
        "RSYLTPGD246N": "R246N"
    }
    m_2 = []
    try:
        mut_list = mut_list.split(",")
        for i in mut_list:
            if i in indel_mapper:
                m_2.append(indel_mapper[i])
            elif ("-" in i) or ("del" in i) or ("ins" in i):
                continue
            else:
                try:
                    int(i[1:-1])
                    m_2.append(i)
                except:
                    continue
    except:
        pass
    return (",".join(m_2))


#### Read in EVEscape scores
evescape_smm = pd.read_csv(EVEscape_path)
evescape_smm["mutations"] = evescape_smm.wt + evescape_smm.i.astype(
    str) + evescape_smm.mut

evescape_smm["evescape_pos"] = evescape_smm["evescape"] - evescape_smm[
    "evescape"].min()
evescape_smm["evescape_sigmoid"] = evescape_smm["evescape"].apply(sigmoid)

# read in GISAID processed data
gisaid = pd.read_csv(gisaid_path)
gisaid["missense mutations"] = gisaid.mutations.apply(clean_up_mutations)
gisaid["number of missense mutations"] = gisaid[
    "missense mutations"].str.count(",") + 1
print("# unique spike sequences >10 times", len(gisaid.query("count >10")))

# Score all sequences
gisaid["EVEscape score_sigmoid"] = gisaid.apply(
    lambda x: sum(evescape_smm[evescape_smm.mutations.isin(x[
        "missense mutations"].split(","))]["evescape_sigmoid"].values),
    axis=1)

gisaid["EVEscape score_pos"] = gisaid.apply(
    lambda x: sum(evescape_smm[evescape_smm.mutations.isin(x[
        "missense mutations"].split(","))]["evescape_pos"].values),
    axis=1)

gisaid.sort_values("EVEscape score_pos", ascending=False, inplace=True)

gisaid.to_csv("../results/summaries_with_gisaid/all_seqs_EVEscape_scores.csv")
