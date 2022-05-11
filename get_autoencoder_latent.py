
######
### PyTorch implementation of an autoencoder trained to approximate Depmap profiles.
### 
### Modified from code written by AFAgarap at https://gist.github.com/AFAgarap/4f8a8d8edf352271fa06d85ba0361f26 
######

# Imports packages
import os
import argparse
import math
import time
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data_utils
import numpy as np
import pandas as pd
from torchsummary import summary
from datetime import timedelta

# Imports model classes
import ae_classes as AE

# Makes parser for user input
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='autoencoder_interface.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file",
                        default=argparse.SUPPRESS,
                        type=str,
                        required=True,
                        help="Path to training data [required]")
    parser.add_argument("-o", "--output_folder",
                        default=argparse.SUPPRESS,
                        type=str,
                        required=True,
                        help="Path to output folder [required]")
    parser.add_argument("-l", "--latent_dim",
                        default=3,
                        type=int,
                        help="Number of latent dimensions for autoencoder [optional]")
    return parser

# Gets latent space from trained autoencoder and writes to file
def get_ae_latent(model, train, n_features, latent_dim, output_folder):

    # Gets latent space corresponding to each profile
    latent_space = np.zeros((train.shape[0], latent_dim))
    fake_profiles = np.zeros((train.shape[0], n_features))
    all_cor = []
    ordered_train_loader = torch.utils.data.DataLoader(torch.from_numpy(train), batch_size=1, shuffle=False)
    i = 0
    print("Train rows: " + str(train.shape[0]))
    print("Train cols: " + str(train.shape[1]))
    for batch in ordered_train_loader:

        # Reshapes profile to size [N, 1, n_features] 
        batch = batch.view(-1, n_features).float()
        latent_space[i,:] = model.get_latent(batch, return_final = True).cpu().detach().numpy()

        # Gets correlation between reconstructed data and real Depmap data
        orig = train[i,:]
        reconstructed = model(batch).cpu().detach().numpy().reshape(orig.shape)
        print(orig.shape)
        print(reconstructed.shape)
        quit()
        cor = np.corrcoef(orig, reconstructed, rowvar = False)[1,0]
        all_cor.append(cor)
        fake_profiles[i,:] = model(batch).cpu().detach().numpy()

        # Increments counter
        i += 1

    # Writes latent space and correlations to file
    np.savetxt(os.path.join(output_folder, "autoencoder_latent_mat.tsv"), latent_space, delimiter="\t")
    np.savetxt(os.path.join(output_folder, "autoencoder_latent_cor.txt"), all_cor)
    np.savetxt(os.path.join(output_folder, "autoencoder_fake_profiles.tsv"), fake_profiles, delimiter="\t")

# Main script
def latent_main(input_file, output_folder, latent_dim):

    # Sets precision for printing
    torch.set_printoptions(precision = 20)

    # Reads in and formats data
    data = pd.read_csv(input_file, sep = "\t", index_col = 0)
    nrow = data.shape[0]
    ncol = data.shape[1]
    train = np.asarray(data)
    n_features = train.shape[1]

    # Sets seed for reproducibility
    seed = 5001

    # Loads model
    model_file = os.path.join(output_folder, "autoencoder_model.pt")
    model = AE.ConvAE(input_shape=n_features, latent_dim=latent_dim)
    model.load_state_dict(torch.load(model_file))
    model.eval()

    # Gets latent space from trained model and writes to file
    get_ae_latent(model, train, n_features, latent_dim, output_folder)
    

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    args = vars(args)
    latent_main(args['input_file'], args['output_folder'], args['latent_dim'])

