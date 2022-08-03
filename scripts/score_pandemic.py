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

#Downloaded metadata.tsv on 5/23/22 from GISAID (must get GISAID access to download)
##Run generate_strains_table.py on the metadata to get the summary_matrix.csv
gisaid_strains = "../data/gisaid/summary_matrix.csv"

##mode lineage per combination of spike mutations from metadata
gisaid_summary_pango = "../data/gisaid/pango_modes.csv"

##spike mutation occurrence by month (from summary_matrix.csv using get_single_mut_dates.py)
count_by_month = "../data/gisaid/single_mutant_count_by_month.csv"

#cumulative single mutation counts (downloaded 5/11/22 from covidcg)
mutation_freq = "../data/gisaid/covidcg_mutation_frequencies.csv"

#consensus mutations for pango lineages (downloaded 7/1/22 from covidcg)
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

#merge mutation counts
rbd_one = rbd_one.merge(freqs[["pos", "ref", "alt", "counts"]],
                        how="left",
                        left_on=["i", "wt", "mut"],
                        right_on=["pos", "ref", "alt"])

#Save RBD table with GISAID counts
rbd_one = rbd_one.drop(columns=["trans", "min_dist", "pos", "ref", "alt"])
rbd_one.to_csv(
    "../results/summaries_with_gisaid/rbd_dist_one_scores_gisaid.csv")

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
for voc_name in ["B.1.1.7", "B.1.351", "B.1.617.2", "P.1", "BA.1"]:
    voc_muts = pango[pango.name == voc_name][pango.alt != "-"][
        pango.ref.str.len() == 1][pango.alt.str.len() ==
                                  1].mutation_name.str.split(
                                      ":", expand=True)[1].values
    spike_one[voc_name] = [x in voc_muts for x in spike_one.mutation]

##Save Spike table with GISAID counts and VOC info
spike_one = spike_one.drop(columns=["trans", "min_dist", "pos", "ref", "alt"])
spike_one.to_csv(
    "../results/summaries_with_gisaid/spike_dist_one_scores_gisaid.csv",
    index=False)

##############################################
#Get data for evescape quantile per month
##############################################

#Read in counts
counts = pd.read_csv(count_by_month)

#Convert to cumulative counts
counts[counts.columns[1:]] = counts[counts.columns[1:]].cumsum(axis=1)

#Reshape
counts = pd.melt(counts,
                 id_vars=["mutation"],
                 var_name="date",
                 value_name="counts")

#Add on EVEscape scores
counts = counts.merge(spike_one[["mutation", "evescape"]],
                      how="left",
                      on="mutation")

counts["evescape_quant"] = [stats.percentileofscore(spike_one.evescape,x)/100 for x in counts.evescape]

#Get evescape quantile bins
counts["evescape_quantile"] = pd.cut(
    counts.evescape_quant, 4, labels=["bottom", "top 75%", "top 50%", "top 25%"])

#number of unique mutations in quantile with counts > 100 before date
freq_plot = []
f = counts[counts.counts > 100]
for index, group in f.groupby(['evescape_quantile', "date"]):
    t, d = index
    freq_plot.append({
        'evescape_quantile':
        t,
        "date":
        d,
        'counts':
        f.loc[(f['date'] <= d)
              & (f['evescape_quantile'].eq(t))]['mutation'].nunique()
    })
freq_plot = pd.DataFrame(freq_plot)

freq_plot.to_csv(
    "../results/summaries_with_gisaid/evescape_quantile_per_month.csv",
    index=False)

##############################################
#Get data for strains to use (combinations occurred >500x)
##############################################
strains = pd.read_csv(gisaid_strains)
strains = strains[strains["count"] > 500]

#Clean up sequence dates
strains["dates"] = [[date for date in x.split(",") if "-" in date]
                    for x in strains.collection_dates.values]
strains["dates"] = [[date if len(date) > 7 else f"{date}-28" for date in x]
                    for x in strains.dates.values]
strains["dates"] = [[dt.strptime(date, '%Y-%m-%d') for date in x]
                    for x in strains.dates.values]

#Take first date (use quantile to avoid year-only dates)
strains["min_date"] = [
    pd.to_datetime(pd.Series(x)).quantile(0.01) for x in strains.dates
]
strains["min_date"] = [
    datetime.date(x.year, x.month, x.day) for x in strains.min_date
]

#Only look at sequences with mutations
strains = strains[~strains.num_mutations.isna()]
strains["mutations"] = strains["mutations"].astype(str)

##ignore insertions and deletions
strains["muts"] = [[
    mut for mut in x.split(",") if "del" not in mut and "ins" not in mut
    and "stop" not in mut and "-" not in mut
] for x in strains.mutations.values]
strains["muts"] = [[mut for mut in x if mut[1:-1].isnumeric()]
                   for x in strains.muts.values]

strains["num_mutations"] = strains.muts.str.len()
strains = strains[strains.num_mutations != 0]

#Get the most common lineage assigned to a unique sequence
modes = pd.read_csv(gisaid_summary_pango)
strains = strains.merge(
    modes.rename(columns={"Pango lineage": "mode_lineage"}),
    left_on="mutations",
    right_on="spike_mutations",
    how="left")

#Pick one GISAID accession ID to assign to sequence
strains["accession_id"] = strains.accession_ids.str.split(",").str.get(0)

##############################################
#Get random combinations of mutations at each depth for zscores
#From set of mutations seen >1000x in GISAID
##############################################

#find mutation counts seen >1000x
freqs["mutation"] = freqs.ref + freqs.pos.astype(str) + freqs.alt
freqs = freqs[freqs.counts > 1000]
freqs = freqs[freqs.alt != "-"][freqs.ref.str.len() == 1][freqs.alt.str.len()
                                                          == 1]
gisaid_muts_1000 = freqs.mutation.values

num = 10000

zscores = pd.DataFrame({"depth": [], "mean": [], "std": []})

for depth in set(strains.num_mutations.values):
    if depth == 1:
        #Use single mutant scores for depth 1
        zscore_seq_scores = spike_one[spike_one.mutation.isin(
            gisaid_muts_1000)]["evescape"].values
        mean = np.mean(zscore_seq_scores)
        std = np.std(zscore_seq_scores)
        zscores.loc[len(zscores.index)] = [depth, mean, std]
    else:
        zscore_seq_scores = []
        seqs_depth = []
        while len(seqs_depth) < num:
            #Make random sets of mutations at depth
            mut_list = list(
                np.random.choice(gisaid_muts_1000, depth, replace=False))
            mut_list.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
            #Check that mutations are at different sites and unique
            if len(set([x[1:-1] for x in mut_list
                        ])) == depth and ":".join(mut_list) not in seqs_depth:
                seqs_depth.append(":".join(mut_list))
                #Sum evescape score for all mutants and append to score list
                score = sum(spike_one[spike_one.mutation.isin(mut_list)]
                            ["evescape"].values)
                zscore_seq_scores.append(score)
        #Get distribution parameters
        mean = np.mean(zscore_seq_scores)
        std = np.std(zscore_seq_scores)
        zscores.loc[len(zscores.index)] = [depth, mean, std]


##############################################
#Get strain-level scores
##############################################
def z_score(score, depth):
    mean = zscores[zscores.depth == depth]["mean"].values[0]
    std = zscores[zscores.depth == depth]["std"].values[0]
    return st.norm(mean, std).cdf(score)


strains["evescape"] = strains.apply(lambda x: sum(spike_one[
    spike_one.mutation.isin(x["muts"])]["evescape"].values),
                                    axis=1)

strains["zscore"] = strains.apply(
    lambda x: z_score(x['evescape'], x['num_mutations']), axis=1)

strains["frequent"] = [x > 50000 for x in strains["count"]]

strains["classification"] = [
    "VOC" if (mode_lineage in voc or mode_lineage in vbm)
    and strains["frequent"].values[index] else ""
    for index, mode_lineage in enumerate(strains.mode_lineage)
]

strains["classification"] = [
    "D614G" if num_mut == 1 and strains["frequent"].values[index] else
    strains["classification"].values[index]
    for index, num_mut in enumerate(strains.num_mutations)
]
strains["classification"] = [
    "other" if strains["classification"].values[index] == "" else
    strains["classification"].values[index]
    for index, mode_lineage in enumerate(strains.mode_lineage)
]

strains = strains[[
    "accession_id", "min_date", "mode_lineage", "num_mutations", "count",
    "evescape", "zscore", "classification", "frequent"
]]
strains.to_csv("../results/summaries_with_gisaid/strain_scores.csv",
               index=False)
