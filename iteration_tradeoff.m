% the script tests two approaches to regularized Janssen iterations on a
% simple case of declipping a short audio signal
%
% it performs the reconstruction for several combinations of the number of
% iterations of the Douglas-Rachford sub-algorithm and number of the
% Janssen iterations, and it evaluates the reconstruction using SDR
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

clear
clc
close all

addpath('utils')
filename = 'iteration_tradeoff_01';

%% set the parameters
p       = 512;   % order of the AR process
theta   = 0.20;  % threshold of the hard clipping
maxit   = 1000;  % iterations of the Janssen algorithm
N       = 8192;  % length of the signal
fi      = 0;     % length of the fade-in / fade-out
gammaC  = 0.1;   % parameter for the DR algorithm (regularized coefficients estimation)
gammaS  = 10;    % parameter for the DR algorithm (signal estimation)
DRmaxit = [5 10 25 50 100 250 500 1000 2500 5000];
lambda  = [0 0.0001 0.001 0.01 0.1 1];

%% initialize the data fields
SDRone = NaN(length(lambda), maxit, length(DRmaxit));
SDRtwo = NaN(length(lambda), maxit, length(DRmaxit));
OBJone = NaN(length(lambda), maxit, length(DRmaxit));
OBJtwo = NaN(length(lambda), maxit, length(DRmaxit));

% if the experiment has been stopped, load the data
if isfile(['results/', filename, '.mat'])
    load(['results/', filename, '.mat'])
end

%% load signal, normalize and crop it
signame = 'violin';
[signal, fs] = audioread(['signals/', signame, '.wav']);
if size(signal, 2) > 1
    signal = mean(signal, 2);
end
cosinus = cos(linspace(0, pi/2, fi)').^2;
signal = signal(fs+1:fs+N);
signal = signal/max(abs(signal));
signal(1:fi) = signal(1:fi).*(1-cosinus);
signal(end-fi+1:end) = signal(end-fi+1:end).*cosinus;

%% create the degraded versions and the associated projection operators
masks.R = abs(signal) < theta;
masks.U = signal >= theta;
masks.L = signal <= -theta;
masks.theta = theta;
degraded = max(-theta, min(signal, theta));

for i = 1:length(DRmaxit)
    
    fprintf('DRmaxit = %d (%d of %d)\n', DRmaxit(i), i, length(DRmaxit))
    
    for j = 1:length(lambda)
        
        fprintf('  lambda = %f (%d of %d)\n', lambda(j), j, length(lambda))

        % if the experiment has been stopped, skip what has been already done
        if ~isnan(SDRtwo(j, end, i))
            continue
        end

        %% declip using our regularization
        [declipped, OBJone(j, 1:maxit, i)] = janssen('declipping', degraded, ...
            masks, lambda(j), p, maxit, 'saveall', true, 'DRmaxit', DRmaxit(i), ...
            'verbose', true);

        %% declip using the heuristic regularization
        [declippedGLP, OBJtwo(j, 1:maxit, i)] = janssen('glp', degraded, ...
            masks, lambda(j), p, maxit, 'saveall', true, 'DRmaxit', DRmaxit(i), ...
            'verbose', true);

        %% evaluate and save SDR and objective
        % we use consistent Janssen, thus SNR can be computed on clipped
        % samples only
        for k = 1:maxit
            SDRone(j, k, i) = snr(signal(~masks.R), signal(~masks.R)-declipped(~masks.R, k));
            SDRtwo(j, k, i) = snr(signal(~masks.R), signal(~masks.R)-declippedGLP(~masks.R, k));
        end

        clear declipped declippedGLP

        save(['results/', filename, '.mat'], 'p', 'N', 'SDRone', 'SDRtwo', 'OBJone', 'OBJtwo', ...
            'lambda', 'maxit', 'DRmaxit', 'gammaC', 'gammaS', 'theta')
    end
end