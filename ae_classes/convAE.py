
######
### PyTorch implementation of an autoencoder trained to approximate Depmap profiles.
### 
### Inspired by code written by AFAgarap at https://gist.github.com/AFAgarap/4f8a8d8edf352271fa06d85ba0361f26 
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

# Convolutional autoencoder class
class ConvAE(nn.Module):
    def __init__(self, input_shape, latent_dim, kernel_size = 3):
        super().__init__()

        # Sets precision for printing
        torch.set_printoptions(precision = 10)

        # Sets variables
        self.input_shape = input_shape
        self.latent_dim = latent_dim
        self.kernel_size = kernel_size

        # Calculates sizes for conversion between convolutional and linear layers
        size1 = self.input_shape - (self.kernel_size - 1) + 2
        size2 = math.ceil(size1 / 2) - (self.kernel_size - 1) + 2
        flattened_size = math.ceil(size2 / 2) * 20

        ### Encoder
        # Reshape goes here...
        self.encoder_conv1 = nn.Conv1d(1, 10, kernel_size, padding=1)
        self.encoder_pool1 = nn.MaxPool1d(2, return_indices=True, padding=1)
        self.encoder_conv2 = nn.Conv1d(10, 20, kernel_size, padding=1)
        self.encoder_pool2 = nn.MaxPool1d(2, return_indices=True, padding=1)
        self.encoder_flatten = nn.Flatten()
        self.encoder_output = nn.Linear(flattened_size, latent_dim)

        ### Decoder
        self.decoder_input = nn.Linear(latent_dim, flattened_size, bias = False)
        # Reshape goes here...
        self.decoder_unpool1 = nn.MaxUnpool1d(2, padding=1)
        self.decoder_conv1 = nn.ConvTranspose1d(20, 10, kernel_size, padding=1)
        self.decoder_unpool2 = nn.MaxUnpool1d(2, padding=1)
        self.decoder_conv2 = nn.ConvTranspose1d(10, input_shape, kernel_size, padding=1)
        self.decoder_flatten = nn.Flatten()
        self.decoder_output = nn.Linear(input_shape**2, input_shape)

    def forward(self, profile):

        # Encodes data into latent space
        x = profile.view([profile.size(0), 1, -1])
        x = self.encoder_conv1(x)
        pool_size1 = x.size()
        x, pool_ind1 = self.encoder_pool1(x)
        x = torch.relu(x)
        x = self.encoder_conv2(x)
        pool_size2 = x.size()
        x, pool_ind2 = self.encoder_pool2(x)
        x = torch.relu(x)
        x = self.encoder_flatten(x)
        x = self.encoder_output(x)
        x = torch.tanh(x)

        # Decodes data from latent space
        x = self.decoder_input(x)
        x = x.view([profile.size(0), 20, -1])
        x = self.decoder_unpool1(x, pool_ind2, pool_size2)
        x = self.decoder_conv1(x)
        x = torch.relu(x)
        x = self.decoder_unpool2(x, pool_ind1, pool_size1)
        x = self.decoder_conv2(x)
        x = torch.relu(x)
        x = self.decoder_flatten(x)
        x = self.decoder_output(x)
        return x

    # Gets latent space corresponding to a given input
    def get_latent(self, profile, return_final = True):
        x = profile.view([profile.size(0), 1, -1])
        x = self.encoder_conv1(x)
        pool_size1 = x.size()
        x, pool_ind1 = self.encoder_pool1(x)
        x = torch.relu(x)
        x = self.encoder_conv2(x)
        pool_size2 = x.size()
        x, pool_ind2 = self.encoder_pool2(x)
        x = torch.relu(x)
        x = self.encoder_flatten(x)
        x = self.encoder_output(x)
        if return_final: 
            # temp = x
            # x = torch.relu(x)
            # x = self.decoder_input(x)
            # x = x.view([profile.size(0), 20, -1])
            # x = self.decoder_unpool1(x, pool_ind2, pool_size2)
            # x = self.decoder_conv1(x)
            # x = torch.relu(x)
            # x = self.decoder_unpool2(x, pool_ind1, pool_size1)
            # print(x)
            # x = self.decoder_conv2(x)
            # x = torch.relu(x)
            # x = self.decoder_flatten(x)
            # x = self.decoder_output(x)
            # return temp
            return torch.tanh(x)
        else:
            return x