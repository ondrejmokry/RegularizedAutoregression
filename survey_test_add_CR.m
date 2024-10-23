clear
clc
close all

toolboxpath = 'survey toolbox';
addpath(genpath(toolboxpath))
addpath('utils')

% name of the results file
filename = 'survey_test';

% name of the waveform directory
directory = ['results/', filename];

% methods used
algos = {'inpainting (c)', ...
         'declipping (c, s)', ...
         'glp (c)'};
method = {'inpainting', 'declipping', 'glp'};

%% load data
load(['results/', filename, '.mat'])

%% process
% for each (incosistent) method, load audio, compute CR and update metrics
for i = 1:length(audio_files)
    for j = 1:length(input_SDRs)

        fprintf(['Audio file: ''', audio_files{i} , '.wav''\n']);

        % load audio files
        S = load([directory, '/', audio_files{i}(1:3), '_', num2str(input_SDRs(j), '%02d'), '.mat']);

        % generate clipping masks
        [data, fs] = audioread([toolboxpath, '/Sounds/' audio_files{i} '.wav']);
        [data_clipped, masks, theta, trueSDR, percentage] = clip_sdr(data, input_SDRs(j));
        masks.U = masks.Mh;
        masks.L = masks.Ml;
        masks.R = masks.Mr;
        
        for m = 1:length(lambdaC)
            for n = 1:length(lambdaS)
                for algo = 1:length(algos)

                    % choose only relevant algorithms
                    if algo ~= 2 || lambdaS(n) == Inf
                        continue
                    end

                    fprintf('Algorithm: %s, lambda = [%s, %s]\n', algos{algo}, num2str(lambdaC(m)), num2str(lambdaS(n)))


                    % load the waveform
                    switch algo
                        case 1, declipped = S.(['inp_', num2str(10*m)]); % actually not needed
                        case 2, declipped = S.(['dec_', num2str(10*m + n)]);
                        case 3, declipped = S.(['glp_', num2str(10*m)]); % actually not needed
                    end

                    %% CR
                    tic
                    restored = crossfaded_replace_reliable(declipped, data_clipped, masks.R);
                    reconstructed(algo).times(i, j, m, n) = reconstructed(algo).times(i, j, m, n) + toc; %#ok<*SAGROW>
                    
                    %% SDR
                    reconstructed(algo).SDRs(i, j, m, n) = sdr(data, restored);
                    fprintf('  SDR of the clipped signal: %.2f dB\n', clipped.SDRs(i, j))
                    fprintf('  SDR of the declipped signal: %.2f dB\n', sdr(data, declipped))
                    fprintf('  SDR of the CR signal: %.2f dB\n\n', reconstructed(algo).SDRs(i, j, m, n))

                    %% PEMO-Q ODG
                    % evaluate the reconstructed signal
                    [~, ~, reconstructed(algo).PEMOQs(i, j, m, n), ~] = audioqual_silent(data, restored, fs);

                    %% PEAQ ODG
                    % save the reference signal as wav
                    data_48 = resample(data, 48000, fs);
                    audiowrite('data_48.wav', data_48, 48000);

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
            save(['results/', filename, '_CR.mat'], ...
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
        end
    end
end

%% clean up
delete data_48.wav restored_48.wav