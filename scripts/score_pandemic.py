import pandas as pd
import numpy as np
from scipy import stats
import scipy.stats as st
import json
from seq_utils import *
from datetime import datetime as dt
import datetime

##############################################
# General Paths
##############################################

#Downloaded metadata.tsv on 10/24 from GISAID (must get GISAID access to download)
##Run generate_strains_table.py on the metadata to get the summary_matrix.csv
##spike mutation occurrence by month (from summary_matrix.csv using get_single_mut_dates.py)
count_by_month = "../data/gisaid/single_mutant_count_by_month.csv"

#cumulative single mutation counts (downloaded 10/31/22 from covidcg)
mutation_freq = "../data/gisaid/covidcg_mutation_frequencies.csv"

#consensus mutations for pango lineages (downloaded 10/26/22 from covidcg)
consensus_muts = "../data/gisaid/covidcg_consensus_mutations.json"

#evescape scores for full spike
spike = "../results/summaries_with_scores/full_spike_evescape.csv"

#evescape scores for rbd
rbd = "../results/summaries_with_scores/spike_rbd_evescape.csv"

##############################################
# Merge Tables
##############################################

#Make minimum codon distance dictionary
AA_Codon = syn_cdn_dict(AminoAcid, Codon_AA)
min_dist_dict = create_min_dist_dictionary(AminoAcid, AA_Codon)

##SPIKE
spike = pd.read_csv(spike)

#Compute minimum codon distance between wt and mut for all muts
spike["trans"] = spike.wt + spike.mut
spike["min_dist"] = [min(min_dist_dict[x], 2) for x in spike.trans]

#Consider only mutations a distance of one away from wuhan
spike_one = spike[spike.min_dist == 1]
spike_one["mutation"] = spike_one.wt + spike_one.i.astype(str) + spike_one.mut

#merge mutation counts
freqs = pd.read_csv(mutation_freq)
spike_one = spike_one.merge(freqs[["pos", "ref", "alt", "counts"]],
                            how="left",
                            left_on=["i", "wt", "mut"],
                            right_on=["pos", "ref", "alt"])

##RBD
rbd = pd.read_csv(rbd)

#Compute minimum codon distance between wt and mut for all muts
rbd["trans"] = rbd.wt + rbd.mut
rbd["min_dist"] = [min(min_dist_dict[x], 2) for x in rbd.trans]

#Consider only mutations a distance of one away from wuhan
rbd_one = rbd[rbd.min_dist == 1]
rbd_one["mutation"] = rbd_one.wt + rbd_one.i.astype(str) + rbd_one.mut

#merge mutation counts
rbd_one = rbd_one.merge(freqs[["pos", "ref", "alt", "counts"]],
                        how="left",
                        left_on=["i", "wt", "mut"],
                        right_on=["pos", "ref", "alt"])


##############################################
#Find voc set of mutations from pango lineages
##############################################

#Read in pango lineage json
with open(consensus_muts) as f:
    pango = json.load(f)
pango = pd.DataFrame(pango)
#Select mutations that appear in >10% of spike lineage
pango = pango[pango.protein == "S"][pango.name != "Unassigned"]
pango = pango[pango.fraction > 0.1]
#Expand out mutation lists
pango["mutation"] = pango.mutation_name.str.split(":", expand=True)[1]
#Drop insertions and deletions
pango = pango[pango.alt != "-"][pango.ref.str.len() == 1][pango.alt.str.len()
                                                          == 1]
#Make list of VOC pango lineages
voc = list(set(pango[pango.name.str.contains("AY")].name.values))
voc.extend([
    "B.1.1.7", "B.1.351", "B.1.617.2", "B.1.427", "B.1.429", "B.1.525",
    "B.1.526", "B.1.617.1", "B.1.617.3", "B.1.621", "B.1.621.1", "P.2",
    "B.1.1.529", "XA", "XB", "XC"
])
voc.extend(set(pango[pango.name.str.contains("P.1")].name.values))
voc.extend(set(pango[pango.name.str.contains("B.1.351")].name.values))
voc.extend(set(pango[pango.name.str.contains("BA.1")].name.values))
voc.extend(set(pango[pango.name.str.contains("BA.2")].name.values))
voc.extend(set(pango[pango.name.str.contains("BA.3")].name.values))
voc.extend(set(pango[pango.name.str.contains("BA.4")].name.values))
voc.extend(set(pango[pango.name.str.contains("BA.5")].name.values))
voc.extend(set(pango[pango.name.str.contains("X")].name.values))

voc = set(voc)

#Get list of VOC mutations
voc_mut_list = list(set(pango[pango.name.isin(voc)].mutation.values))

#Label VOC mutations in spike table
spike_one["voc_mut"] = [x in voc_mut_list for x in spike_one.mutation]

##Add columns to indicate if mutation is in selected list of VOCs
for voc_name in ["B.1.1.7", "B.1.351", "B.1.617.2", "P.1", "BA.1", "BA.2","BA.4", "BA.2.12.1","BA.2.75"]:
    voc_muts = pango[pango.name == voc_name][pango.alt != "-"][
        pango.ref.str.len() == 1][pango.alt.str.len() ==
                                  1].mutation_name.str.split(
                                      ":", expand=True)[1].values
    spike_one[voc_name] = [x in voc_muts for x in spike_one.mutation]

    
##################################################################################
#Get first month where each mutation reaches a total of 100 observations in GISAID
##################################################################################

counts = pd.read_csv(count_by_month)

#Convert to cumulative counts
counts[counts.columns[1:]] = counts[counts.columns[1:]].cumsum(axis=1)

#reshape
counts = pd.melt(counts,
                 id_vars=["mutation"],
                 var_name="date",
                 value_name="counts")

#subset to counts>100
counts = counts[counts.counts > 100]

#get first month
counts = counts.sort_values(by="date")
counts = counts.drop_duplicates(subset="mutation")
counts["first_seen"] = counts.date

#merge dates
rbd_one = rbd_one.merge(counts[["mutation","first_seen"]],
                        how="left",
                        on="mutation")
spike_one = spike_one.merge(counts[["mutation","first_seen"]],
                        how="left",
                        on="mutation")


##############################################
#Save Files
##############################################

#Save RBD table with GISAID counts
rbd_one = rbd_one.drop(columns=["trans", "min_dist", "pos", "ref", "alt"])
rbd_one.to_csv(
    "../results/summaries_with_gisaid/rbd_dist_one_scores_gisaid.csv")

##Save Spike table with GISAID counts and VOC info
spike_one = spike_one.drop(columns=["trans", "min_dist", "pos", "ref", "alt"])
spike_one.to_csv(
    "../results/summaries_with_gisaid/spike_dist_one_scores_gisaid.csv",
    index=False)



