clc;
close all;
clear;

addpath('core/');

%% load timeseries for a single subject

% preprocessed and parcellated resting-state data
% 400 x 1200 timeseries for a single subject
% and the network affiliation vector

data = load('./data/hcp_1subj.mat');

v0 = data.v0;               % n x t timeseries
v0 = zscore(v0,[],2);       % zscore data
m0 = data.m0;               % network affiliation vector

%% PCA constraints

k = 2;  % number of components to be constrained

v1 = gradient_sampler(v0, k); % ssample data

% get PCA for empirical and model timeseries
pc0 =  pca(v0', 'NumComponents', k+1);
pc1 =  pca(v1', 'NumComponents', k+1);

figure

for jj = 1:3
    subplot(1, 3, jj);
    plot(pc0(:,jj), pc1(:,jj), '.');
    title(['PC ' num2str(jj)]);
    xlabel('empirical loadings');
    ylabel('model loadings');
    axis equal;
    axis tight;
end
