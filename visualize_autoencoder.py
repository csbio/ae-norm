
######
### Trains autoencoder on Depmap data.
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

# Imports torchviz for visualizing model architecture
from torchviz import make_dot

# Makes parser for user input
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='visualize_autoencoder.py',
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

# Visualizes a single batch
def vizualize_ae(model, train, n_features):
    ordered_train_loader = torch.utils.data.DataLoader(torch.from_numpy(train), batch_size=1, shuffle=False)
    first_batch = next(iter(ordered_train_loader))
    first_batch = first_batch.view(-1, n_features).float()
    make_dot(model(first_batch), params=dict(model.named_parameters()), show_attrs=True, show_saved=True).render("attached", format="pdf")

# Main script
def visualize_main(input_file, output_folder, latent_dim):

    # Sets precision for printing
    torch.set_printoptions(precision = 20)

    # Reads in and formats data
    data = pd.read_csv(input_file, sep = "\t", index_col = 0)
    nrow = data.shape[0]
    ncol = data.shape[1]
    train = np.asarray(data)
    n_features = train.shape[1]
    print("Rows: " + str(nrow))
    print("Cols: " + str(ncol))

    # Loads model
    model_file = os.path.join(output_folder, "autoencoder_model.pt")
    model = AE.ConvAE(input_shape=n_features, latent_dim=latent_dim)
    model.load_state_dict(torch.load(model_file))
    model.eval()

    # Gets latent space from trained model and writes to file
    vizualize_ae(model, train, n_features)
    

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    args = vars(args)
    visualize_main(args['input_file'], args['output_folder'], args['latent_dim'])


