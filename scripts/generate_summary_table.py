import datetime
import getopt
import json
import math
import os
import sys
import pandas as pd
'''
Print the usage message.
'''


def usage():
    print((
        'Usage: python generate_summary_table.py [-h] [--from_date FROM_DATE] '
        '[--to_date TO_DATE] PATH_TO_METADATA_FILE\n'
        'Generate a summary table from a given GISAID metadata file.\n\n'
        'Arguments\n'
        '   PATH_TO_METADATA                        The path to a metadata file. Must be a .tsv. \n\n'
        'Options\n'
        '   --force                                 Force overwrite of output_dir.\n'
        '   -f FROM_DATE, --from_date=FROM_DATE     The lower bound date. Comparisons will be\n'
        '                                           made from this date. Must be YYYY-MM-DD.\n'
        '                                           If not provided, treats all metadata as new.\n'
        '   -t TO_DATE, --to_date=TO_DATE           The upper bound date. Comparisons will be\n'
        '                                           made to this date. Must be YYYY-MM-DD.\n'
        '                                           Defaults to today.\n'
        '   -o DIR, --output_dir=DIR                The directory where output will be stored.\n'
        '                                           Will be created if it doesn\'t exist.\n'
        '                                           Defaults to "./output".'))


'''
From a list of all aa substitutions, extract and sort those prefixed with Spike_
because those are the ones in the spike sequence.
'''


def extract_spike_mutations(mutations):
    # Remove the parentheses and turn mutations into a list
    mutations = mutations[1:-1].split(',')

    spike_mutations = []
    spike_indicator = 'Spike_'
    for mutation in mutations:
        if mutation.startswith(spike_indicator):
            # The mutations aren't sorted, so pull out the index and insert
            # mutations as tuples to be sorted later
            mutation_residue_index = ''.join(c for c in mutation
                                             if c.isdigit())
            # Sometimes it just says "Spike_" with nothing else. Go figure.
            if mutation_residue_index:
                spike_mutations.append((int(mutation_residue_index),
                                        mutation.replace(spike_indicator, '')))
    # Sort based on index, then remove those indices
    sorted_spike_muations = [t[1] for t in sorted(spike_mutations)]

    return ','.join(sorted_spike_muations)


'''
Process the metadata file and return a summary matrix and a subset
of that summary matrix of just new sequences (i.e. sequences that only appear
after from_date)
'''


def process_gisaid_metadata(filepath, from_date, to_date):
    summary_matrix = {}
    nans_seen = {}
    non_conforming_dates_seen = {}
    entries_processed = 0
    entries_since_date_processed = 0

    delimiter = '\t'
    chunksize = 10**5

    # Read in the CSV in chunks. Specify certain datatypes because they are
    # ambiguous when reading in the metadata file.
    with pd.read_csv(filepath,
                     chunksize=chunksize,
                     sep='\t',
                     dtype={
                         'Additional location information': str,
                         'Is reference?': str,
                         'Is complete?': str
                     }) as reader:
        chunks_read = 0
        for chunk in reader:
            print(f'Processing chunk number {chunks_read}')
            chunks_read += 1
            for i in range(len(chunk)):
                entries_processed += 1
                row = chunk.iloc[i]
                if isinstance(row['AA Substitutions'], float) and math.isnan(
                        row['AA Substitutions']):
                    nans_seen['AA Substitutions'] = nans_seen.get(
                        'AA Substitutions', 0) + 1
                    continue
                spike_mutations = extract_spike_mutations(
                    row['AA Substitutions'])
                if not spike_mutations in summary_matrix:
                    summary_matrix[spike_mutations] = {
                        'mutations': spike_mutations,
                        'num_mutations': len(spike_mutations.split(',')),
                        'count': 0,
                        'count_since_date': 0,
                        'accession_ids': [],
                        'collection_dates': [],
                        'pango_lineages': [],
                        'pangolin_versions': [],
                        'nearest_VOCs': [],
                    }
                summary_row = summary_matrix[spike_mutations]

                try:
                    date = datetime.datetime.strptime(row['Collection date'],
                                                      "%Y-%m-%d")
                except ValueError:
                    try:
                        date = datetime.datetime.strptime(
                            row['Collection date'], "%Y-%m")
                        non_conforming_dates_seen[row[
                            'Collection date']] = non_conforming_dates_seen.get(
                                row['Collection date'], 0) + 1
                    except ValueError:
                        try:
                            date = datetime.datetime.strptime(
                                row['Collection date'], "%Y")
                            non_conforming_dates_seen[row[
                                'Collection date']] = non_conforming_dates_seen.get(
                                    row['Collection date'], 0) + 1
                        except ValueError:
                            raise ValueError(
                                f'Unable to convert {row["Collection date"]} to any date format.'
                            )
                            sys.exit(2)

                # If this entry's date is after the period we care about, skip it.
                if date > to_date:
                    continue
                summary_row['count'] += 1

                if date >= from_date:
                    entries_since_date_processed += 1
                    summary_row['count_since_date'] += 1

                summary_row['collection_dates'].append(row['Collection date'])

                # To deal with NaNs, always make sure the values are strings
                if isinstance(row['Accession ID'], str):
                    summary_row['accession_ids'].append(row['Accession ID'])
                elif isinstance(row['Accession ID'], float) and math.isnan(
                        row['Accession ID']):
                    nans_seen['Accession ID'] = nans_seen.get(
                        'Accession ID', 0) + 1

                if (isinstance(
                        row['Pango lineage'],
                        str)) and (row['Pango lineage']
                                   not in summary_row['pango_lineages']):
                    summary_row['pango_lineages'].append(row['Pango lineage'])
                elif isinstance(row['Pango lineage'], float) and math.isnan(
                        row['Pango lineage']):
                    nans_seen['Pango lineage'] = nans_seen.get(
                        'Pango lineage', 0) + 1

                if (isinstance(
                        row['Pangolin version'],
                        str)) and (row['Pangolin version']
                                   not in summary_row['pangolin_versions']):
                    summary_row['pangolin_versions'].append(
                        row['Pangolin version'])
                elif isinstance(row['Pangolin version'], float) and math.isnan(
                        row['Pangolin version']):
                    nans_seen['Pangolin version'] = nans_seen.get(
                        'Pango lineage', 0) + 1

                if ((isinstance(row['Variant'], str)) and
                    (row['Variant'] not in summary_row['nearest_VOCs'])):
                    summary_row['nearest_VOCs'].append(row['Variant'])
                elif isinstance(row['Variant'], float) and math.isnan(
                        row['Variant']):
                    nans_seen['Variant'] = nans_seen.get('Variant', 0) + 1

    new_sequences = {}
    new_sequences_count = 0
    total_sequences_count = 0
    for k in summary_matrix:
        total_sequences_count += 1
        summary_row = summary_matrix[k]
        summary_row['accession_ids'] = ','.join(
            sorted(summary_row['accession_ids']))
        summary_row['pango_lineages'] = ','.join(
            sorted(summary_row['pango_lineages']))
        summary_row['pangolin_versions'] = ','.join(
            sorted(summary_row['pangolin_versions']))
        summary_row['nearest_VOCs'] = ','.join(
            sorted(summary_row['nearest_VOCs']))
        summary_row['proportion_since_date'] = (
            f'{(summary_row["count_since_date"]/entries_since_date_processed):.4f}'
        )

        # Because of non-standard data entry, the 'first_seen' column should
        # be populated in the following priority:
        # 1. The earliest date with the format YYYY-MM-DD
        # 2. If there are no dates in YYYY-MM-DD format, the earliest date in
        #    YYYY-MM format.
        # 3. If there are no dates in YYYY-MM format, the earliest YYYY date.
        y_m_seen = False
        y_seen = False
        collection_dates_sorted = sorted(summary_row['collection_dates'])
        for date in collection_dates_sorted:
            try:
                date = datetime.datetime.strptime(date, "%Y-%m-%d")
                summary_row['first_seen'] = date
                break
            except ValueError:
                if not y_m_seen:
                    try:
                        date = datetime.datetime.strptime(date, "%Y-%m")
                        summary_row['first_seen'] = date
                    except ValueError:
                        if not y_seen:
                            date = datetime.datetime.strptime(date, "%Y")
                            summary_row['first_seen'] = date
        summary_row['collection_dates'] = ','.join(collection_dates_sorted)

        if not isinstance(summary_row['count'], int):
            summary_row['count'] = int(summary_row['count'])
        if not isinstance(summary_row['count_since_date'], int):
            summary_row['count_since_date'] = int(
                summary_row['count_since_date'])
        if summary_row['count'] == summary_row['count_since_date']:
            new_sequences[k] = summary_row
            new_sequences_count += 1

    run_summary = (
        'Run Summary:\n'
        '-------------------\n'
        f'Total number of entries processed: {entries_processed}\n'
        f'Total number of unique sequences seen: {total_sequences_count}\n'
        f'Number of new sequences seen between {from_date.date()} and {to_date.date()}: {new_sequences_count}\n'
        f'Number of NaNs seen: {len(nans_seen)}\n'
        f'Number of non-conforming dates seen: {sum([non_conforming_dates_seen[k] for k in non_conforming_dates_seen])}\n'
        'Dict of NaNs seen:\n'
        f'{json.dumps(nans_seen, sort_keys=True, indent=2)}\n'
        'Dict of non-conforming dates seen:\n'
        f'{json.dumps(non_conforming_dates_seen, sort_keys=True, indent=2)}\n')
    print(run_summary)
    return pd.DataFrame.from_dict(summary_matrix).T, pd.DataFrame.from_dict(
        new_sequences).T, run_summary


if __name__ == "__main__":
    # Process options and arguments
    from_date, to_date = '', ''
    force = False
    output_dir = './output'
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], 'hf:t:o:',
            ['help', 'force', 'from_date=', 'to_date=', 'output_dir='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('--force'):
            force = True
        elif opt in ('-f', '--from_date'):
            from_date = arg
        elif opt in ('-t', '--to_date'):
            to_date = arg
        elif opt in ('-o', '--output_dir'):
            output_dir = arg.rstrip('/')

    if ((os.path.exists(output_dir)) and (not force)):
        sys.exit(f'Output directory {output_dir} already exists. '
                 'Please either use the --force flag or specify a different '
                 'directory using -o/--output_dir.')

    if from_date:
        try:
            from_date = datetime.datetime.strptime(from_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                "Error: Incorrect format for FROM_DATE. Must be YYYY-MM-DD.")
            sys.exit(2)
    else:
        from_date = datetime.datetime.min

    if to_date:
        try:
            to_date = datetime.datetime.strptime(to_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                "Error: Incorrect format for TO_DATE. Must be YYYY-MM-DD.")
    else:
        to_date = datetime.datetime.today()

    if len(args) != 1:
        sys.exit('Please provide the path to a GISAID metadata file!')
    filepath = args[0]
    if not os.path.exists(filepath):
        sys.exit(f'File at {filepath} doesn\'t exist.')
    elif not filepath.endswith('.tsv'):
        sys.exit(f'Please provide a *.tsv file.')

    summary_matrix, new_sequences, run_summary = process_gisaid_metadata(
        filepath, from_date, to_date)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    summary_matrix.to_csv(f'{output_dir}/summary_matrix.csv', index=False)
    new_sequences.to_csv(f'{output_dir}/new_sequences.csv', index=False)
    with open(f'{output_dir}/run_summary.txt', 'w') as f:
        f.write(run_summary)
