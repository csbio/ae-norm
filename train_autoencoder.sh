#!/bin/bash

# Important parameters
num_epochs=1
epochs_str="epochs_${num_epochs}"
input_file=input/depmap_q2_2020_nona_mean_rst_clp_mms.tsv
output_folder="output/autoencoder/q2_2020_rst_clp_mms/$epochs_str/cnn_latent_"

# Static paths
olf_file=input/static/olfactory_receptors.csv
essentials_file=input/static/depmap_common_essentials_2019_q3.txt
corum_file=input/static/corum_human.tsv

# Runs pipeline for different latent space values
for latent_dim in 1 2 3 4 5 10 25
do
	latent_folder="$output_folder$latent_dim"
	pr_folder="${latent_folder}/pr"
	mkdir -p "$latent_folder"
	python3 train_autoencoder.py -i "$input_file" -o "$latent_folder" -e $num_epochs -l "$latent_dim"
	python3 get_autoencoder_latent.py -i "$input_file" -o "$latent_folder" -l "$latent_dim"
	Rscript pr_prep_autoencoder.R -i "$input_file" -f "$olf_file" -e "$essentials_file" -c "$corum_file" -o "$latent_folder" -n $latent_dim
	Rscript run_flex_autoencoder.R -o "$latent_folder" -m
	Rscript analyze_complexes.R -o "$pr_folder" -e "$essentials_file" -c "$corum_file"
done
