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

%% Network models

nc0 = net_corr(v0, m0);      % empirical network correlation

% Intra-network correlations
v1 = network_sampler(v0, m0, 'intra');
nc1 = net_corr(v1, m0);      % model network correlation

% All-network correlations
v2 = network_sampler(v0, m0, 'all');
nc2 = net_corr(v2, m0);      % model network correlation

figure
subplot(121)
plot(nc0, nc1, '.');
xlabel('empirical network correlations');
ylabel('model network correlations');
title('intra-network model')
axis equal;
axis square;

subplot(122)
plot(nc0, nc2, '.');
xlabel('empirical network correlations');
ylabel('model network correlations');
title('all-network model')
axis equal;
axis square;

