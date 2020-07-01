%%%% Demo of BLOTLESS dictionary learning
%%%% Demostration of BLOTLESS dictionary learning algorithm.
clear all
close all
clc
addpath(genpath('utilities'))
%% Generate a problem
m = 64;                        % dictionary atom dimension
l = 128;                       % number of dictionary atoms, l>m means over-complete dictionary
k = 5;                         % Column sparsity
n = 600;                       % Number of samples
D = randn(m,l);                
D = column_normalize(D,2);     % Column-normalized Dictionary

%% Sparse Subspace
X = gen_coef(l,n,k) ;          % generate random sparse coefficients
Y = D*X;
% Y = awgn(Y,20,'measured');   % Noisy or not 

%% Initial D
D_ini = randn(m,l);            % Initial Dictionary
D_ini = column_normalize(D_ini,2);
iter_num = 100;                % Overall iterations

%% Parameters
params.Tdata = k;              % sparisity level
params.iternum = iter_num;
params.initialDictionary = D_ini;
params.groudtruthDict = D;
params.verbose = 'yes';

%% BLOTLESS Dictionary Learning
[D_hat,Rec_err] = Blotless_orig_func(Y,params);

%% Draw figure
figure;plot(1:iter_num,Rec_err);xlabel('Iteration Number');ylabel('Dictionary Recovery Error')



