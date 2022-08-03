import pandas as pd
import numpy as np
import scipy
from scipy import stats

##############################################
#Escape Thresholds
##############################################

bloomthresh = 0.57
xiethresh = 0.9

##############################################
#General paths
##############################################

#evescape scores for rbd by site
rbd_sites = "../results/summaries_with_scores/spike_rbd_evescape_sites.csv"

#escape dms for voc sera (downloaded from https://github.com/jbloomlab/SARS-CoV-2-RBD_Delta/blob/main/results/supp_data/aggregate_raw_data.csv)
voc_escape_sera = '../data/experiments/bloom_rbd_escape/aggregate_raw_data_strains.csv'

#evescape scores for rbd
rbd_scores = "../results/summaries_with_scores/spike_rbd_evescape.csv"

#rbd with all experiments
rbd_experiments = "../results/summaries/rbd_experiments_and_scores.csv"

#antibody properties from https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP_Vir_mAbs/blob/main/data/antibody_annotations_input.csv
ab_properties = "../data/antibody_properties/antibody_annotations_input.csv"

#antibody metadata processed from escape calculator data file
ab_metadata = "../data/antibody_properties/rbd_antibody_metadata.csv"

##############################################
#Wuhan and VOC sera escape comparison
##############################################

rbd_sites = pd.read_csv(rbd_sites)

#calculate sera escape for sera from different strains
for lib in ["Delta", "Beta", "Wuhan-Hu-1"]:
    new_rbd_dat = pd.read_csv(voc_escape_sera)

    #remove monoclonal antibodies from table
    new_rbd_dat = new_rbd_dat[~new_rbd_dat["class"].str.contains("class")]

    #use only primary delta/delta reinfection for delta
    if lib == "Delta":
        new_rbd_dat = new_rbd_dat[new_rbd_dat["class"].str.contains("Delta")]

    new_rbd_dat = new_rbd_dat[new_rbd_dat.library == lib]
    new_rbd_dat["mutation"] = new_rbd_dat.wildtype + new_rbd_dat.site.astype(
        str) + new_rbd_dat.mutation
    dat = new_rbd_dat.pivot(index='mutation',
                            columns='condition',
                            values=['mut_escape'])

    #mean site total sera
    dat["site"] = dat.index.str[1:-1]
    dat["site"] = dat["site"].astype(int)
    dat = dat.groupby("site").sum()
    rbd_sites = rbd_sites.merge(
        dat.mean(axis=1).to_frame(), left_on="i", right_on="site",
        how="left").rename(columns={0: f"mean_escape_sera_{lib}"})

#get max across strains
rbd_sites["max_sera"] = rbd_sites[[
    "mean_escape_sera_Wuhan-Hu-1", "mean_escape_sera_Delta",
    "mean_escape_sera_Beta"
]].max(axis=1)

#get max strain for plot hue
#wuhan if above threshold otherwise max if other voc above threshold
rbd_sites["max_sera_strain"] = rbd_sites[[
    "mean_escape_sera_Wuhan-Hu-1", "mean_escape_sera_Delta",
    "mean_escape_sera_Beta"
]].idxmax(axis=1)

rbd_sites["max_sera_strain"] = rbd_sites["max_sera_strain"].mask(
    (rbd_sites["mean_escape_sera_Wuhan-Hu-1"] > bloomthresh),
    "mean_escape_sera_Wuhan-Hu-1")

rbd_sites["max_sera_strain"] = rbd_sites["max_sera_strain"].mask(
    (rbd_sites["max_sera"] <= bloomthresh), "No Escape")

rbd_sites.to_csv(
    "../results/summaries_with_added_dms/spike_rbd_evescape_voc_sera.csv",
    index=False)

##############################################
#Assign each mutation to antibody class
##############################################

##read in mutation level df with experiments
rbd_experiments = pd.read_csv(rbd_experiments)

#add binarized escape to experiment results
rbd = pd.read_csv(rbd_scores)
rbd_experiments = rbd_experiments.merge(rbd[[
    "wt", "mut", "i", 'is_escape_experiment_bloom', "is_escape_experiment_xie"
]],
                                        on=["wt", "mut", "i"],
                                        how="left")

rbd_experiments["mutation"] = rbd_experiments.wt + rbd_experiments.i.astype(
    str) + rbd_experiments.mut
ab_metadata = pd.read_csv(ab_metadata)

#get clases for each antibody
classes_dict = dict(zip(ab_metadata.condition, ab_metadata.condition_subtype))

#Get subset of escape experiments that are antibodies in Bloom data
ab_conditions = set(
    ab_metadata[ab_metadata.condition_type == "antibody"].condition.values)

rbd_experiments_conditions = [
    x for x in rbd_experiments.columns.values
    if "escape_" in x and "_Bloom" in x and x in ab_conditions
]

#For each mutant, get list of Bloom antibodies that it escapes
rbd_experiments['escapes_bloom'] = rbd_experiments[
    rbd_experiments_conditions].apply(lambda x: list(rbd_experiments[
        rbd_experiments_conditions].columns[x > bloomthresh]),
                                      axis=1)

rbd_experiments['escapes_bloom'] = [[] if x is np.NaN else x
                                    for x in rbd_experiments['escapes_bloom']]

#Get subset of escape experiments that are antibodies in Xie data
rbd_experiments_conditions = [
    x for x in rbd_experiments.columns.values
    if "escape_" in x and "_Xie" in x and x in ab_conditions
]

#For each mutant, get list of Xie antibodies that it escapes
rbd_experiments['escapes_xie'] = rbd_experiments[
    rbd_experiments_conditions].apply(lambda x: list(rbd_experiments[
        rbd_experiments_conditions].columns[x > xiethresh]),
                                      axis=1)

rbd_experiments['escapes_xie'] = [[] if x is np.NaN else x
                                  for x in rbd_experiments['escapes_xie']]

#get total list of antibodies escaped by mutation
rbd_experiments["escapes"] = (rbd_experiments["escapes_bloom"] +
                              rbd_experiments["escapes_xie"])

#return the most common antibody class of all antibodies escaped by mutation
rbd_experiments["class"] = [[classes_dict[condition] for condition in set(x)]
                            for x in rbd_experiments.escapes]
rbd_experiments["class"] = rbd_experiments["class"].apply(
    lambda x: scipy.stats.mode(np.array(x))[0])
rbd_experiments["class"] = [
    "" if len(x) == 0 else x[0] for x in rbd_experiments["class"]
]

##############################################
#Antibody properties
##############################################

#read in antibody properties
ab_properties = pd.read_csv(ab_properties)
ab_properties = ab_properties[~ab_properties.value.isna()]

#select antibodies with breadth and neutralization data
bloom_abs = set([
    x for x in rbd_experiments.columns.values
    if "escape_" in x and "_Bloom" in x
])
breadth_ab = ab_properties[ab_properties.metric.isin(["breadth"
                                                      ])].antibody.values
neut_ab = ab_properties[ab_properties.metric.isin(["SARS2_IC50", "Zost_IC50"
                                                   ])].antibody.values
rbd_ab_conditions_with_neut = [
    x for x in bloom_abs
    if x.split("_")[1:-1] in breadth_ab and x.split("_")[1:-1] in neut_ab
]

#get max escape for antibodies with breadth and neutralization
rbd_experiments['max_escape_with_neut'] = rbd_experiments[
    rbd_ab_conditions_with_neut].max(axis=1)
rbd_experiments['max_escape_exp_with_neut'] = rbd_experiments[
    rbd_ab_conditions_with_neut].idxmax(axis=1)
rbd_experiments[
    "max_escape_exp_with_neut"] = rbd_experiments.max_escape_exp_with_neut.str.split(
        "_", expand=True)[1]

#add onto main rbd df with evescape scores
rbd = rbd.merge(rbd_experiments[[
    "wt",
    "mut",
    "i",
    "class",
    'max_escape_with_neut',
    "max_escape_exp_with_neut",
]],
                on=["wt", "mut", "i"],
                how="left")

#add antibody properties
rbd["breadth"] = [
    ab_properties[ab_properties.antibody == ab][ab_properties.metric ==
                                                "breadth"].value.values[0]
    for ab in rbd["max_escape_exp_with_neut"]
]
rbd["neutralization"] = [
    ab_properties[ab_properties.antibody == ab][ab_properties.metric.isin(
        ["SARS2_IC50", "Zost_IC50"])].value.values[0]
    for ab in rbd["max_escape_exp_with_neut"]
]
rbd["neutralization"] = np.log10(rbd["neutralization"])

rbd.to_csv(
    "../results/summaries_with_added_dms/spike_rbd_evescape_antibody_properties.csv",
    index=False)
