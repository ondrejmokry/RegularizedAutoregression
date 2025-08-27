% tested algorithms:
% (1)  inpainting
% (2)  inpainting (c)
% (3)  consistent declipping
% (4)  consistent declipping (c)
% (5)  consistent declipping (s)
% (6)  consistent declipping (c, s)
% (7)  inconsistent declipping
% (8)  inconsistent declipping (c)
% (9)  inconsistent declipping (s)
% (10) inconsistent declipping (c, s)
% (11) glp
% (12) glp (c)
% applied with lambda = 0.01 and lambda = [0.01, 10] for the
% consistent and inconsistent variants, respectively

% tested variants:
% (1) 10 main iterations, 1000 iterations of DR
% (2) 10 main iterations, progressive iterations of DR
% (3) 5 main iterations, 1000 iterations of DR, linesearch
% (4) 5 main iterations, 1000 iterations of DR, extrapolation of the
%     estimated AR coefficients
% (5) 5 main iterations, 1000 iterations of DR, extrapolation of the
%     estimated signal
% (6) 5 main iterations, 1000 iterations of DR, extrapolation of the
%     estimated AR coefficients and the estimated signal

% tested pw settings:
% p = [512, 2048]
% w = [2048, 4096, 8192]

% tested input SDRs: [ 5 7 10 15 ] dB

% tested signals:
% (1) a08_violin
% (2) a18_bassoon
% (3) a35_glockenspiel
% (4) a42_accordion
% (5) a60_piano_schubert
% from signals/small_set.mat

% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

clear
clc
close all

addpath(genpath('survey toolbox'))
addpath('utils')

filename = 'acceleration_test_rect_02';

%% load audio files
signals = load('signals/small_set.mat');
% signals = load('survey toolbox/Sounds/Sound_database.mat');

%% possible parameter values
audio_files = signals.names(1:2:end);
% audio_files = fieldnames(signals);
% signals.fs = 44100;
% input_SDRs  = [ 5, 7, 10, 15 ];
input_SDRs  = 10;
ps          = [ 512, 2048 ];
ws          = [ 2048, 4096, 8192 ];
variants    = 6;
algos       = {'inpainting', ...
               'inpainting (c)', ...
               'consistent declipping', ...
               'consistent declipping (c)', ...
               'consistent declipping (s)', ...
               'consistent declipping (c, s)', ...
               'inconsistent declipping', ...
               'inconsistent declipping (c)', ...
               'inconsistent declipping (s)', ...
               'inconsistent declipping (c, s)', ...
               'glp', ...
               'glp (c)'};

%% algos and variants settings
% 12 algos
method    = {'inpainting', 'inpainting', ...
             'declipping', 'declipping', 'declipping', 'declipping', 'declipping', 'declipping', 'declipping', 'declipping', ...
             'glp', 'glp'};
coefaccel = [0 1 0 1 0 1 0 1 0 1 0 1];
sigaccel  = [0 0 0 0 1 1 0 0 2 2 0 0];
lambda    = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, ...
            [0.01, 10], [0.01, 10], [0.01, 10], [0.01, 10], ...
            0.01, 0.01};

% 6 variants
maxit      = [10, 10, 5, 5, 5, 5];
DRmaxit    = {1000, round(logspace(2, 3, 10))', 1000, 1000, 1000, 1000};
linesearch = [0 1 0 0 0 0];
coefextra  = [0 0 0 1 0 1];
sigextra   = [0 0 0 0 1 1];

%% initialization of the data fields
% clipped signals
clipped.SDRs   = NaN(length(audio_files), length(input_SDRs));
clipped.PEMOQs = NaN(length(audio_files), length(input_SDRs));
clipped.PEAQs  = NaN(length(audio_files), length(input_SDRs));

% reconstructed signals
reconstructed = struct();
for algo = 1:length(algos)
    reconstructed(algo).algorithm = algos{algo};
    reconstructed(algo).SDRs   = NaN(length(audio_files), length(input_SDRs), length(ps), length(ws), variants);
    reconstructed(algo).PEMOQs = NaN(length(audio_files), length(input_SDRs), length(ps), length(ws), variants);
    reconstructed(algo).PEAQs  = NaN(length(audio_files), length(input_SDRs), length(ps), length(ws), variants);
    reconstructed(algo).times  = NaN(length(audio_files), length(input_SDRs), length(ps), length(ws), variants);
end

% number of combinations
combs = numel(reconstructed(1).SDRs);
comb  = 0;
cinit = datetime("now");

% if the experiment has been stopped, load the data
if isfile(['results/', filename, '.mat'])
    load(['results/', filename, '.mat'])
end
skipped = 0;

for i = 1:length(audio_files)
    for j = 1:length(input_SDRs)
        for k = 1:length(ps)
            for l = 1:length(ws)
                for m = 1:variants
                    
    comb = comb + 1;
    
    % if the experiment has been stopped, skip what has been already done
    if ~isnan(reconstructed(end).SDRs(i, j, k, l, m))
        skipped = skipped + 1;
        continue
    end
    
    % estimate the remaining time
    c = datetime("now");
    d = seconds(c-cinit); % elapsed time in seconds
    d = d/3600; % elapsed time in hours
    fprintf('======================================\n')
    fprintf('Combination %d of %d\n', comb, combs)
    fprintf('Elapsed time: %d hours\n', round(d))
    fprintf('Estimated remaining time: %d hours\n', ...
    	round(d*((combs-skipped)/(comb-skipped)-1)))
    fprintf('======================================\n')
    
    %% load audio-file and clip it
    % load
    data = signals.(audio_files{i});
    
    % clip (using the function in survey toolbox/Tools)
    [data_clipped, masks, theta, trueSDR, percentage] = clip_sdr(data, input_SDRs(j));

    % save masks in the right form for segmentation
    masks.U = masks.Mh;
    masks.L = masks.Ml;
    masks.R = masks.Mr;

    for algo = 1:length(algos)
        %% main algorithm
        fprintf('Algorithm: %s (%d of %d)\n', algos{algo}, algo, length(algos))
        tic
        restored = janssen(...
            method{algo}, data_clipped, masks, lambda{algo}, ps(k), maxit(m), ... % model and main parameters
            'segmentation', true, 'wtype', 'rect', 'w', ws(l), 'a', ws(l)/2, ...  % overlap-add
            'coefaccel', coefaccel(algo), 'sigaccel', sigaccel(algo), ...         % acceleration
            'DRmaxit', DRmaxit{m}, ...                                            % inner iterations
            'linesearch', linesearch(m), ...                                      % liensearch
            'coefextra', coefextra(m), 'sigextra', sigextra(m), ...               % extrapolation
            'verbose', false);                                                    % (not) printing processed segment
        reconstructed(algo).times(i, j, k, l, m) = toc;
        
        %% SDR
        clipped.SDRs(i, j) = sdr(data, data_clipped);
        reconstructed(algo).SDRs(i, j, k, l, m) = sdr(data, restored);
        fprintf('  SDR of the clipped signal: %.2f dB\n', clipped.SDRs(i, j))
        fprintf('  SDR of the reconstructed signal: %.2f dB\n', reconstructed(algo).SDRs(i, j, k, l, m))

        %% PEMO-Q ODG
        % evaluate the clipped signal
        if k == 1 && l == 1 && m == 1
            [~, ~, clipped.PEMOQs(i, j), ~] = audioqual_silent(data, data_clipped, signals.fs);
        end
        
        % evaluate the reconstructed signal
        [~, ~, reconstructed(algo).PEMOQs(i, j, k, l, m), ~] = audioqual_silent(data, restored, signals.fs);

        %% PEAQ ODG
        if k == 1 && l == 1 && m == 1
            % save the reference signal as wav
            data_48 = resample(data, 48000, signals.fs);
            audiowrite('data_48.wav', data_48, 48000);
        
            % save the clipped signal as wav
            clipped_48 = resample(data_clipped, 48000, signals.fs);
            audiowrite('clipped_48.wav', clipped_48, 48000);
            
            % evaluate the clipped signal
            [clipped.PEAQs(i, j), ~] = PQevalAudio_fn('data_48.wav', 'clipped_48.wav', 0, length(clipped_48));
        end
        
        % save the restored signal as wav
        restored_48 = resample(restored, 48000, signals.fs);
        audiowrite('restored_48.wav', restored_48, 48000);
        
        % evaluate the restored signal
        [reconstructed(algo).PEAQs(i, j, k, l, m), ~] = PQevalAudio_fn('data_48.wav', 'restored_48.wav', 0, length(restored_48));
    end
    
    %% save
    d = datetime("now");
    save(['results/', filename, '.mat'], ...
        'algos', ...
        'audio_files', ...
        'clipped', ...
        'd', ...
        'input_SDRs', ...
        'ps', ...
        'reconstructed', ...
        'ws')
    
                end
            end
        end
    end
end

%% clean up
delete data_48.wav clipped_48.wav restored_48.wav