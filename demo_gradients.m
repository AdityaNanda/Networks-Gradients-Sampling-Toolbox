clc;
close all;
clear;

addpath('core/');

%% load timeseries for a single subject

% preprocessed and parcellated resting-state data
% v0, timeseries of a single subject (400 x 1200)
% m0, network affiliation vector     (400 x 1)

data = load('./data/hcp_1subj.mat');

v0 = data.v0;                   % timeseries
v0 = zscore(v0,[],2);           % zscore data
m0 = data.m0;                   % network affiliation

%% PCA constraints

k = 2;  % number of components to be constrained

% Sample data
v1 = gradient_sampler(v0, k);

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
