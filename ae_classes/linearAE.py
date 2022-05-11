
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

# Linear autoencoder class
class LinearAE(nn.Module):
    def __init__(self, input_shape, latent_dim):
        super().__init__()
        self.encoder_hidden = nn.Linear(in_features=input_shape, out_features=latent_dim)
        self.encoder_output = nn.Linear(in_features=latent_dim, out_features=latent_dim)
        self.decoder_hidden = nn.Linear(in_features=latent_dim, out_features=latent_dim)
        self.decoder_output = nn.Linear(in_features=latent_dim, out_features=input_shape)

    def forward(self, features):
        activation = self.encoder_hidden(features)
        activation = torch.relu(activation)
        code = self.encoder_output(activation)
        code = torch.relu(code)
        activation = self.decoder_hidden(code)
        activation = torch.relu(activation)
        activation = self.decoder_output(activation)
        reconstructed = torch.tanh(activation)
        return reconstructed

    def get_latent(self, features):
        activation = self.encoder_hidden(features)
        activation = torch.relu(activation)
        code = self.encoder_output(activation)
        code = torch.relu(code)
        return(code)