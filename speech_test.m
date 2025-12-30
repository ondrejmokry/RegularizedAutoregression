% the script performs the declipping experiment similar to the survey
% 
%    P. Záviška, P. Rajmic, A. Ozerov and L. Rencker, "A Survey and an
%    Extensive Evaluation of Popular Audio Declipping Methods, " in IEEE
%    Journal of Selected Topics in Signal Processing, vol. 15, no. 1, pp.
%    5-24, Jan. 2021, doi: 10.1109/JSTSP.2020.3042071.
%
% for several settings of the two AR-based algorithms inspired by Janssen:
%
%    (1) modification of the original Janssen algorithm for inpainting
%        which includes the AR coefficients regularization
%    (2) modification of the original Janssen algorithm which includes the
%        consistency constraints and AR coefficients regularization
%    (3) the method from the US patent US8126578B2:
%        Clipped-waveform repair in acoustic signals using generalized
%        linear prediction (https://patents.google.com/patent/US8126578)
%        with AR coef regularization (will be referred to as GLP)
%
% the difference is that this script loads speech signals instead of music
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

clear
clc
close all

toolboxpath = 'survey toolbox';
addpath(genpath(toolboxpath))
addpath('utils')

% name of the file to save the results to
filename = 'speech_test';

% name of the directory to save the waveforms to
directory = ['results/', filename];
mkdir(directory)

%% possible parameter values
audio_files = dir('speech');
audio_files = audio_files(3:end-1);

input_SDRs = [ 3, 5, 7, 10 ];
p          = 256;  % AR model order
w          = 2048; % window length
lambdaC    = [ 0, 1e-5, 1e-3 ];
lambdaS    = [ 10, Inf ];
algos      = {'inpainting (c)', ...
              'declipping (c, s)', ...
              'glp (c)'};
method     = {'inpainting', 'declipping', 'glp'};

%% initialization of the data fields
% clipped signals
clipped.SDRs   = NaN(length(audio_files), length(input_SDRs));
clipped.PEMOQs = NaN(length(audio_files), length(input_SDRs));
clipped.PEAQs  = NaN(length(audio_files), length(input_SDRs));

% reconstructed signals
reconstructed = struct();
for algo = 1:length(algos)
    reconstructed(algo).algorithm = algos{algo};
    if algo == 2
        dims = [length(audio_files), length(input_SDRs), length(lambdaC), length(lambdaS)];
    else
        dims = [length(audio_files), length(input_SDRs), length(lambdaC)];
    end
    reconstructed(algo).SDRs   = NaN(dims);
    reconstructed(algo).PEMOQs = NaN(dims);
    reconstructed(algo).PEAQs  = NaN(dims);
    reconstructed(algo).times  = NaN(dims);
end

% number of combinations
combs = numel(reconstructed(1).SDRs)+numel(reconstructed(2).SDRs)+numel(reconstructed(3).SDRs);
comb  = 0;
cinit = datetime("now");

% if the experiment has been stopped, load the data
if isfile(['results/', filename, '.mat'])
    load(['results/', filename, '.mat'])
end
skipped = 0;

for i = 1:length(audio_files)
    for j = 1:length(input_SDRs)

        %% load audio-file and clip it
        % load
        fprintf(['Loading audio ''', audio_files(i).name, '\n']);
        [data, fs] = audioread(['speech/' audio_files(i).name]);

        % clip (using the function in survey toolbox/Tools)
        [data_clipped, masks, theta, trueSDR, percentage] = clip_sdr(data, input_SDRs(j));

        % save masks in the right form for segmentation.m
        masks.U = masks.Mh;
        masks.L = masks.Ml;
        masks.R = masks.Mr;
        
        for m = 1:length(lambdaC)
            for n = 1:length(lambdaS)
                for algo = 1:length(algos)

                    % skip lambdaS iteration for inpainting and glp
                    if algo ~= 2 && n > 1
                        continue
                    end

                    comb = comb + 1;

                    % if the experiment has been stopped, skip what has been already done
                    if ~isnan(reconstructed(algo).SDRs(i, j, m, n))
                        skipped = skipped + 1;
                        continue
                    end

                    % estimate the remaining time
                    c = datetime('now');
                    d = seconds(c-cinit); % elapsed time in seconds
                    d = d/3600; % elapsed time in hours
                    fprintf('Algorithm: %s, lambda = [%s, %s] (combination %d of %d)\n', algos{algo}, num2str(lambdaC(m)), num2str(lambdaS(n)), comb, combs)
                    fprintf('Elapsed time: %d hours\n', round(d))
                    fprintf('Estimated remaining time: %d hours\n', ...
                        round(d*((combs-skipped)/(comb-skipped)-1)))

                    %% main algorithm
                    tic
                    restored = janssen(...
                        method{algo}, data_clipped, masks, [lambdaC(m), lambdaS(n)], p, 5, ... % model and main parameters
                        'segmentation', true, 'wtype', 'rect', 'w', w, 'a', w/2, ...           % overlap-add
                        'coefaccel', true, 'sigaccel', true, ...                               % acceleration
                        'DRmaxit', 1000, ...                                                   % inner iterations
                        'linesearch', false, ...                                               % liensearch
                        'coefextra', false, 'sigextra', true, ...                              % extrapolation
                        'verbose', false);                                                     % (not) printing processed segment
                    reconstructed(algo).times(i, j, m, n) = toc;
                    
                    % save the waveform to the struct (to be saved)
                    switch algo
                        case 1, S.(['inp_', num2str(10*m)]) = restored;
                        case 2, S.(['dec_', num2str(10*m + n)]) = restored;
                        case 3, S.(['glp_', num2str(10*m)]) = restored;
                    end

                    %% SDR
                    clipped.SDRs(i, j) = sdr(data, data_clipped);
                    reconstructed(algo).SDRs(i, j, m, n) = sdr(data, restored);
                    fprintf('  SDR of the clipped signal: %.2f dB\n', clipped.SDRs(i, j))
                    fprintf('  SDR of the reconstructed signal: %.2f dB\n\n', reconstructed(algo).SDRs(i, j, m, n))

                    %% PEMO-Q ODG
                    % evaluate the clipped signal
                    if max([m, n, algo]) == 1
                        [~, ~, clipped.PEMOQs(i, j), ~] = audioqual_silent(data, data_clipped, fs);
                    end

                    % evaluate the reconstructed signal
                    [~, ~, reconstructed(algo).PEMOQs(i, j, m, n), ~] = audioqual_silent(data, restored, fs);

                    %% PEAQ ODG
                    if max([m, n, algo]) == 1
                        % save the reference signal as wav
                        data_48 = resample(data, 48000, fs);
                        audiowrite('data_48.wav', data_48, 48000);

                        % save the clipped signal as wav
                        clipped_48 = resample(data_clipped, 48000, fs);
                        audiowrite('clipped_48.wav', clipped_48, 48000);

                        % evaluate the clipped signal
                        [clipped.PEAQs(i, j), ~] = PQevalAudio_fn('data_48.wav', 'clipped_48.wav', 0, length(clipped_48));
                    end

                    % save the restored signal as wav
                    restored_48 = resample(restored, 48000, fs);
                    audiowrite('restored_48.wav', restored_48, 48000);

                    % evaluate the restored signal
                    [reconstructed(algo).PEAQs(i, j, m, n), ~] = PQevalAudio_fn('data_48.wav', 'restored_48.wav', 0, length(restored_48));

                end
            end

            %% save
            % metrics
            d = datetime('now');
            save(['results/', filename, '.mat'], ...
                'algos', ...
                'audio_files', ...
                'clipped', ...
                'd', ...
                'input_SDRs', ...
                'lambdaC', ...
                'lambdaS', ...
                'p', ...
                'reconstructed', ...
                'w')

            % waveforms
            S.clean = data;
            savename = [directory, '/', audio_files(i).name(1:end-4), '_', num2str(input_SDRs(j), '%02d'), '.mat'];
            if isfile(savename)
                save(savename, '-struct', 'S', '-append')
            else
                save(savename, '-struct', 'S')
            end
            clear S
            
        end
    end
end

%% clean up
delete data_48.wav clipped_48.wav restored_48.wav