# -*- coding: utf-8 -*-
"""
@File   ： GCN_VAE_frequency.py
@Time   ： 2023/1/3 21:30
@Author ： Jia Yiming
"""
import argparse
import numpy as np
import torch
from matplotlib import pyplot as plt
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv
import torch.nn as nn



class GCN_VAE_t(torch.nn.Module):
    def __init__(self, latent_size = 10):
        super().__init__()
        self.latent_size = latent_size
        self.conv1 = GCNConv(50, 50)
        self.encoder_forward1 = nn.Sequential(
            nn.Linear(50, 40),
            nn.LeakyReLU(),
            nn.Linear(40, 30),
            nn.LeakyReLU(),
            nn.Linear(30, 20)
        )

        self.FC = nn.Sequential(
            nn.Linear(10, 10),
            nn.ReLU()
        )
        self.decoder_forward1 = nn.Sequential(
            nn.Linear(10, 20),
            nn.LeakyReLU(),
            nn.Linear(20, 40),
            nn.LeakyReLU(),
            nn.Linear(40, 50),
            nn.Sigmoid()
        )


    def encoder(self, X, encoder_forward):
        out = encoder_forward(X)
        mu = out[:, :self.latent_size]
        log_var = out[:, self.latent_size:]
        return mu, log_var

    def decoder(self, z, decoder_forward):
        mu_prime = decoder_forward(z)
        return mu_prime


    def reparameterization(self, mu, log_var):
        epsilon = torch.randn_like(log_var)
        z = mu + epsilon * torch.sqrt(log_var.exp())
        return z

    def loss(self, X, mu_prime, mu, log_var):
        # reconstruction_loss = F.mse_loss(mu_prime, X, reduction='mean') is wrong!
        reconstruction_loss = torch.mean(torch.square(X - mu_prime).sum(dim=1))

        # latent_loss = torch.mean(0.5 * (log_var.exp() + torch.square(mu) - log_var).sum(dim=1))

        # loss = reconstruction_loss + latent_loss

        loss = reconstruction_loss

        return loss


    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x1 = x[:, 100:]
        x1 = self.conv1(x1, edge_index)
        mu1, log_var1 = self.encoder(x1, self.encoder_forward1)
        z1 = self.reparameterization(mu1, log_var1)
        z = self.FC(z1)
        mu_prime1 = self.decoder(z, self.decoder_forward1)
        r1 = [mu_prime1, mu1, log_var1]

        return r1, z

def train(model, optimizer, data_loader, device, name='GCN_VAE_t'):
    model.train()

    total_loss = 0
    for X in data_loader:
        # batch_size = X.shape[0]
        X = X.to(device)
        model.zero_grad()
        r1, z = model(X)

        loss = model.loss(X.x[:, 100:], r1[0], r1[1], r1[2])

        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss


def main(args):
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    device = torch.device(device)

    batch_size = args.batch_size
    epochs = args.n_epoch
    lr = args.lr
    edge_index = torch.tensor(np.load("graph_index.npy").transpose(), dtype=torch.long)

    x = torch.tensor(np.float32(np.load('res1.npy')))
    data = Data(x=x, edge_index=edge_index)
    # data.edge_index = add_self_loops(data.edge_index)[0]


    data_loader = DataLoader(dataset=[data], batch_size=1, shuffle=False)

    # train VAE
    gcn_vae = GCN_VAE_t().to(device)
    optimizer = torch.optim.AdamW(gcn_vae.parameters(), lr=lr)

    all_loss = []
    all_epoch = np.array(range(epochs))
    print('Start Training VAE...')
    for epoch in range(1, 1 + epochs):
        loss = train(gcn_vae, optimizer, data_loader, device, name='GCN_VAE_f')
        all_loss.append(loss)
        print("Epochs: {epoch}, AvgLoss: {loss:.4f}".format(epoch=epoch, loss=loss))
    print('Training for VAE has been done.')

    np.array(all_loss)
    plt.plot(all_epoch, all_loss)
    plt.show()

    PATH = 'GCN_VAE_f.pth'
    torch.save(gcn_vae.state_dict(), PATH)
    gcn_vae.eval()

    _,  z = gcn_vae(data)

    z = z.cpu().detach().numpy()
    print(z.shape)
    np.save('mu_f', z)


if __name__ == '__main__':
    print("____________________________")
    # parameters setting
    parser = argparse.ArgumentParser()
    parser.add_argument("--lr", type=float, default=0.001, help="learnin rate")
    parser.add_argument("--n_epoch", type=int, default=100, help="number of epoch")
    parser.add_argument("--hidden_dim", type=int, default=20, help="dimention of hidden layers")
    parser.add_argument("--batch_size", type=int, default=32, help="size of batch")
    args = parser.parse_args()
    main(args)
