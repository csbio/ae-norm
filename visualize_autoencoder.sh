#!/bin/bash

# Important parameters
num_epochs=1
epochs_str="epochs_${num_epochs}"
input_file=input/depmap_q2_2020_nona_mean_rst_clp_mms.tsv
output_folder="output/autoencoder/q2_2020_rst_clp_mms/$epochs_str/cnn_latent_"

# Runs pipeline for different latent space values
for latent_dim in 1
do
	latent_folder="$output_folder$latent_dim"
	pr_folder="${latent_folder}/pr"
	mkdir -p "$latent_folder"
	python3 visualize_autoencoder.py -i "$input_file" -o "$latent_folder" -l "$latent_dim"
done
