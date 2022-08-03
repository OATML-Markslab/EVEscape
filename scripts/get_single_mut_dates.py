import pandas as pd
import datetime

##############################################
# General Paths
##############################################

##Downloaded metadata.tsv on 5/23/22 from GISAID (must get GISAID access to download)
##Run generate_strains_table.py on the metadata to get the summary_matrix.csv
gisaid_strains = "../data/gisaid/summary_matrix.csv"

##############################################
#Fucntions to aggregate by month
##############################################


def create_month_list(end_month):
    months = []
    for yr in ["2020", "2021"]:
        for m in range(1, 13):
            if m < 10:
                months.append(yr + "-0" + str(m))
            else:
                months.append(yr + "-" + str(m))

    for yr in ["2022"]:
        for m in range(1, end_month + 2):
            months.append(yr + "-0" + str(m))
    return (months)


def process_dates(row):
    mutation_list = row.replace("[", "").replace("]",
                                                 "").replace("'", "").replace(
                                                     " ", "").split(",")
    mutation_list = [
        i for i in mutation_list if i not in ["2020", "2021", "2022"]
    ]  #remove year only mutations
    return (mutation_list)


def month_counts_dict(mutation_list):
    month_counts = dict([(key, 0) for key in months])
    for i in mutation_list:
        for m in range(len(months) - 1):
            if months[m] <= i <= months[m + 1]:
                month_counts[months[m]] += 1
    return (month_counts)


##############################################
#Get single mutation counts per month
##############################################

#make months list
months = create_month_list(5)

#read in metadata summary file
gisaid = pd.read_csv(gisaid_strains)

#get set of single mutations
all_muts = list(
    gisaid.mutations.dropna().str.split(","))  #all mutations in GISAID
concat_muts = [item for sublist in all_muts for item in sublist]  #flatten list
concat_muts = set(concat_muts)

#filter mutations
df = pd.DataFrame({"mutations": list(concat_muts)})
df = df[~df["mutations"].str.contains("ins|del|stop|B|X|Z")]

## collect dates
dates = []
for i, mut in enumerate(list(df.mutations)):
    mut_subset = gisaid[gisaid['mutations'].str.contains(mut, na=False)]
    all_dates = list(mut_subset.collection_dates.dropna().str.split(
        ","))  #all dates with mutation
    concat_dates = [item for sublist in all_dates
                    for item in sublist]  #flatten list
    dates.append(concat_dates)
    if i % 100 == 0:
        print(i, flush=True)

#save to df
df = pd.DataFrame({"mutation": list(df.mutations), "dates": dates})
df["site"] = [int(m[1:-1]) for m in df.mutation]

#aggregate months
df_100 = df[df.dates.str.len() > 99]  # only keep sequences seen more than 100x

df_100["dates"] = df_100["dates"].astype(str)
df_100["dates"] = df_100.dates.apply(process_dates)
df_100["month_counts"] = df_100.dates.apply(month_counts_dict)

df_fin = pd.DataFrame.from_records(list(df_100["month_counts"]))
df_fin.insert(0, "mutation", list(df_100["mutation"]))

#save to csv
df_fin.to_csv("../data/gisaid/single_mutant_count_by_month.csv", index=False)
