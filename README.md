# EVEscape

This is the official code repository for the paper "Learning from pre-pandemic data to forecast viral antibody escape" (). This paper is a joint collaboration between the Marks lab (https://www.deboramarkslab.com/) and the OATML group (https://oatml.cs.ox.ac.uk/).

##Overview
EVEscape is a model that computes the predicted likelihood of a given viral protein variant to induce immune escape from antibodies. For each protein, EVEscape predicts escape from data sources available pre-pandemic: sequence likelihood predictions from broader viral evolution, antibody accessibility information from protein structures, and changes in binding interaction propensity from residue chemical properties.     

##Usage
Computing EVEscape scores consists of three components:
1. Fitness: use scores from EVE, an unsupervised generative model of mutation effect from broader evolutionary sequences   
2. Accessibility: calculate WCN from PDB structures of relevant conformations of the viral protein of interest
3. Dissimilarity: calculate difference in charge and hydrophobicity between the mutant residue and the wildtype 

The components are then standardized and fed into a temperature scaled logistic function, and we take the the log transform of the product of the 3 terms to obtain final EVEscape scores. 

We also provide EVEscape scores for all single mutation variants of SARS-CoV-2 Spike and aggregate strain-level predictions for all PANGO lineages in our paper, and EVEscape rankings of newly occurring GISAID ranking and visualization of likely future mutations will be available at evescape.org. 

##Example scripts
The scripts folder contains python scripts to calculate EVEscape scores for all single mutations and aggregate deep mutational scanning data for SARS-CoV-2 RBD, Flu HA, and HIV Env from [./data]. 
Specifically this includes the following two scripts:
 - [scripts/process_protein_data.py] calculates the three EVEscape components 
 - [scripts/evescape_scores.py] creates the final evescape scores and outputs scores and processed DMS data in [./results/summaries_with_scores]

##Data requirements
The data required to obtain EVEscape scores is one or multiple PDB files and EVE scores (see [https://github.com/OATML-Markslab/EVE] for how to generate) and a fasta file of the wildtype sequence for the viral protein of interest. 

##License
This project is available under the MIT license. 

##Reference
If you use this code, please cite the following paper:

Nicole N. Thadani*, Sarah Gurev*, Pascal Notin*, Noor Youssef, Nathan J. Rollins, Chris Sander, Yarin Gal, Debora S. Marks. Learning from pre-pandemic data to forecast viral antibody escape. BioRxiv. 2022. (* equal contribution)

Links:
 - Pre-print: 
 - Website: evescape.org
