
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
    parser.add_argument("-e", "--epochs",
                        default=10,
                        type=int,
                        help="Number of training epochs [optional]")
    parser.add_argument("-l", "--latent_dim",
                        default=3,
                        type=int,
                        help="Number of latent dimensions for autoencoder [optional]")
    return parser

# Trains Depmap on CPUs with Adam optimizer and MSE loss
def train_ae(train, n_features, epochs, latent_dim, output_folder):

    # Prints start time
    start_time = time.time()
    print(output_folder)
    with open(os.path.join(output_folder, "time.txt"), "w") as f:
        print("Start time: " + str(start_time), file = f)

    # Converts numpy array to tensor
    train_loader = torch.utils.data.DataLoader(torch.from_numpy(train), batch_size=32, shuffle=True)

    # Sets up optimizer
    device = torch.device("cpu")
    model = AE.ConvAE(input_shape=n_features, latent_dim=latent_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.MSELoss()

    # Trains for the given number of epochs
    all_loss = []
    for epoch in range(epochs):
        loss = 0
        for batch in train_loader:
            
            # Reshapes mini-batch to size [N, 1, n_features] 
            batch = batch.view(-1, n_features).to(device).float()
            
            # Resets gradients
            optimizer.zero_grad()
            
            # Computes reconstructions, loss, and gradients
            outputs = model(batch)
            train_loss = criterion(outputs, batch)
            train_loss.backward()
            
            # Updates parameters and appends loss
            optimizer.step()
            loss += train_loss.item()
        
        # Computes and prints the average epoch training loss
        loss = loss / len(train_loader)
        all_loss.append(loss)
        print("epoch : {}/{}, loss = {:.6f}".format(epoch + 1, epochs, loss))

    # Writes epoch loss to file
    np.savetxt(os.path.join(output_folder, "autoencoder_mse_loss.txt"), all_loss)

    # Prints end time and time taken
    end_time = time.time()
    with open(os.path.join(output_folder, "time.txt"), "w") as f:
        print("End time: " + str(end_time), file = f)
        print("Time taken: " + str(timedelta(seconds = end_time - start_time)), file = f)

    # Returns trained model
    return model

# Main script
def train(input_file, output_folder, epochs, latent_dim):

    # Sets precision for printing
    torch.set_printoptions(precision = 100)

    # Reads in and formats data
    data = pd.read_csv(input_file, sep = "\t", index_col = 0)
    nrow = data.shape[0]
    ncol = data.shape[1]
    train = np.asarray(data)
    n_features = train.shape[1]

    # Sets seed for reproducibility
    seed = 5001

    # Trains autoencoder and saves model
    model_file = os.path.join(output_folder, "autoencoder_model.pt")
    model = train_ae(train, n_features, epochs, latent_dim, output_folder)
    torch.save(model.state_dict(), model_file)
    

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    args = vars(args)
    train(args['input_file'], args['output_folder'], args['epochs'], args['latent_dim'])

