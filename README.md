# EVEscape

This is the official code repository for the paper ["Learning from pre-pandemic data to forecast viral escape"](https://www.nature.com/articles/s41586-023-06617-0). This paper is a joint collaboration between the [Marks Lab](https://www.deboramarkslab.com/) and the [OATML group](https://oatml.cs.ox.ac.uk/).

## Overview
EVEscape is a model that computes the predicted likelihood of a given viral protein variant to induce immune escape from antibodies. For each protein, EVEscape predicts escape from data sources available pre-pandemic: sequence likelihood predictions from broader viral evolution, antibody accessibility information from protein structures, and changes in binding interaction propensity from residue chemical properties.   

## Usage
Computing EVEscape scores consists of three components:
1. Fitness: use scores from EVE, an unsupervised generative model of mutation effect from broader evolutionary sequences   
2. Accessibility: calculate WCN from PDB structures of relevant conformations of the viral protein of interest
3. Dissimilarity: calculate difference in charge and hydrophobicity between the mutant residue and the wildtype 

The components are then standardized and fed into a temperature scaled logistic function, and we take the the log transform of the product of the 3 terms to obtain final EVEscape scores. 

We also provide EVEscape scores for all single mutation variants of SARS-CoV-2 Spike and aggregate strain-level predictions for all GISAID strains, and EVEscape rankings of newly occurring GISAID strains and visualization of likely future mutations will be available at evescape.org. 

## Scripts
The scripts folder contains python scripts to calculate EVEscape scores for all single mutations and aggregate available deep mutational scanning data for SARS-CoV-2 RBD, Flu HA, HIV Env, Lassa glycoprotein, and Nipah fusion and glycoproteins from [data](/data). 
Specifically this includes the following two scripts:
 - [process_protein_data.py](scripts/process_protein_data.py) calculates the three EVEscape components 
 - [evescape_scores.py](scripts/evescape_scores.py) creates the final evescape scores and outputs scores and processed DMS data in [summaries_with_scores](./results/summaries_with_scores)
 
 The scripts folder also contains a python script [score_pandemic_strains.py](scripts/score_pandemic_strains.py) to calculate EVEscape scores for all strains in GISAID. The output strain scores (~150MB unzipped) can be downloaded as follows:
 ```
curl -o strain_scores_20230318.zip https://marks.hms.harvard.edu/evescape/strain_scores_20230318.zip
unzip strain_scores_20230318.zip
rm strain_scores_20230318.zip
```
The workflow of the scripts to create the data tables in [results](./results) needed for the main figures of the EVEscape paper is available in [evescape_summary.pdf](./evescape_summary.pdf). Additional data tables are available in the paper supplement. 

## Data requirements
The training data required to obtain EVEscape scores is one or multiple PDB files of the viral antigen, EVE scores (see next subsection) and a fasta file of the wildtype sequence for the viral protein of interest. 

To download the RBD escape data used as validation data in this project (~120MB unzipped):
```
curl -o escape_dms_data_20220109.zip https://marks.hms.harvard.edu/evescape/escape_dms_data_20220109.zip
unzip escape_dms_data_20220109.zip
rm escape_dms_data_20220109.zip
```
(originally downloaded from [SARS2_RBD_Ab_escape_maps](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps))

## Generating EVE scores
We leverage the original [EVE codebase](https://github.com/OATML-Markslab/EVE) to compute the evolutionary indices used in EVEscape.

### Model training
The MSAs used to train the EVE models used in this project can be found in the supplemental material of the paper (Data S1). 

We modify the Bayesian VAE [training script](https://github.com/OATML-Markslab/EVE/blob/master/train_VAE.py) to support the following hyperparameter choices in the [MSA_processing](https://github.com/OATML-Markslab/EVE/blob/master/utils/data_utils.py) call:
- sequence re-weighting in MSA (theta): we choose a value of 0.01 that is better suited to viruses (Hopf et al., Riesselman et al.)
- fragment filtering (threshold_sequence_frac_gaps): we keep sequences in the MSA that align to at least 50% of the target sequence.
- position filtering (threshold_focus_cols_frac_gaps): we keep columns with at least 70% coverage, except for SARS-CoV-2 Spike for which we lower the required value to 30% in order to maximally cover experimental positions and significant pandemic sites.

We train 5 independent models with different random seeds.

### Model scoring
For the 5 independently-trained models, we compute [evolutionary indices](https://github.com/OATML-Markslab/EVE/blob/master/compute_evol_indices.py) sampling 20k times from the approximate posterior distribution (ie., num_samples_compute_evol_indices=20000). We then average the resulting scores across the 5 models to obtain the final EVE scores used in EVEscape.

### Model checkpoints
We provide open access to the EVE models we trained for the various viruses discussed in the paper. These model checkpoints were obtained following the training procedure described above. To download checkpoints for a viral protein of interest, please adapt the following example with the relevant filename (filenames are listed in the table underneath).
```
curl -o EVE_checkpoints_I4EPC4.zip https://marks.hms.harvard.edu/evescape/EVE_checkpoints_I4EPC4.zip
unzip EVE_checkpoints_I4EPC4.zip
rm EVE_checkpoints_I4EPC4.zip
```
| Uniprot ID     | Organism          | Protein name         | Filename          |
| :---------------- | :---------------- | :---------------- | :---------------- | 
| I4EPC4            |Influenza A virus   | Hemagglutinin   | EVE_checkpoints_I4EPC4.zip |
| P0DTC2 (pre2020 sequences only) |SARS-CoV-2         | Spike glycoprotein  | EVE_checkpoints_P0DTC2_full_pre2020.zip |
| P0DTC2     |SARS-CoV-2         | Spike glycoprotein    | EVE_checkpoints_P0DTC2_full.zip |
| Q2N0S5      | HIV-1       | Envelope glycoprotein gp160    | EVE_checkpoints_Q2N0S5.zip |
| FUS_NIPAV    | Nipah virus   | Fusion glycoprotein F0  | EVE_checkpoints_FUS_NIPAV.zip |
| GLYC_LASSJ   | Lassa virus  | Pre-glycoprotein polyprotein GP complex | EVE_checkpoints_GLYC_LASSJ.zip |
| GLYCP_NIPAV   | Nipah virus  |  Glycoprotein G   | EVE_checkpoints_GLYCP_NIPAV.zip |


## Software requirements
The entire codebase is written in python. The corresponding environment may be created via conda and the provided [requirements.txt](./requirements.txt) file as follows:
```
conda config --add channels conda-forge
conda create --name evescape_env --file requirements.txt
conda activate evescape_env
```
The environment installs in minutes.

## Runtime
After collecting the training data, generating EVEscape scores for all single mutations runs in minutes. Strain scoring of all GISAID strains runs in 2 hours on 64G of memory. 

## License
This project is available under the MIT license. 

## Reference
If you use this code, please cite the following paper:

Nicole N. Thadani*, Sarah Gurev*, Pascal Notin*, Noor Youssef, Nathan J. Rollins, Daniel Ritter, Chris Sander, Yarin Gal, Debora S. Marks. Learning from pre-pandemic data to forecast viral escape. _Nature_. 2023. 

(* equal contribution)

Links:
 - Publication: https://www.nature.com/articles/s41586-023-06617-0
 - Website: https://www.evescape.org/

See new work using EVEscape to design infectious Spike proteins that forecast future neutralizing antibody escape on [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.10.08.561389v1).
