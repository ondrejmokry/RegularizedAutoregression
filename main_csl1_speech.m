%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%           (P(W))CSL1           %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%          (all sounds)          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the declipping methods according to B. Defraene et al: 
% Declipping of audio signals using perceptual compressed sensing
% 
% Pavel Záviška, Brno University of Technology, 2020
% adapted by Ondøej Mokrý for the speech dataset, 2025
% 
% using toolbox LTFAT
ltfatstart

addpath('survey toolbox/Tools')
addpath(genpath('survey toolbox/Evaluation_algorithms'))
addpath('survey toolbox/Methods/PCSL1')

close all
clear variables

%% Database
sounds = dir('speech');
sounds = sounds(3:end-1);

algorithms = {'csl1', 'pcsl1', 'pwcsl1'};
input_SDRs = [ 3, 5, 7, 10 ];

% initialization of matrices for SDR and dSDR values, and computational time
SDR = NaN(length(sounds), length(input_SDRs));
dSDR = NaN(length(sounds), length(input_SDRs));
TIME = NaN(length(sounds), length(input_SDRs));

% initialization of counter
cnt = 0;

for sound = 1:length(sounds)
    for clip_idx = 1:length(input_SDRs)
        for algorithm = 1
            
            cnt = cnt+1;
            %% input file settings
            % load audio-file
            [data, fs] = audioread(['speech/' sounds(sound).name]);
            
            % signal length
            param.Ls = length(data);
            
            % sampling frequency
            param.fs = fs; % Hz

            
            %% settings
            param.inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold
            
            % window parameters
            param.w = 2048;       % window length
            param.a = param.w/4;  % window shift
            param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html
            
            % DFT parameters
            param.F = frame('dft');
            param.F.redundancy = 2;  % redundancy of the DFT transform
            param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
            param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);
            
            % algorithm
            param.algorithm = algorithms{algorithm}; % algorithm to compute declipping, options: 'aspade', 'sspade', 'sspade_new'

            % paramsolver parameters
            paramsolver = settings_csl1(param.algorithm);
            
            paramsolver.store_dsdr = 0; 
            paramsolver.store_obj = 0;
            
            paramsolver.verbose = false;
                                              
            %% clipping
            [data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
            
            
            %% Main algorithm
            tic;
            
            [data_rec, sdr_iter, obj_iter] = csl1_segmentation(data_clipped, param, paramsolver, data);
            
            time = toc;
            
            
            %% Time & SDR evaluation
            
            % Time
            TIME(sound, clip_idx) = time;

            % Save restored file
            if input_SDRs(clip_idx) < 10
                eval([sounds(sound).name(1:end-4) '_rec_' algorithms{algorithm} '_0' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            else
                eval([sounds(sound).name(1:end-4) '_rec_' algorithms{algorithm} '_' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            end

            % SDR
            sdr_clip = sdr(data, data_clipped);
            sdr_rec = sdr(data, data_rec);

            SDR(sound, clip_idx) = sdr_rec;
            dSDR(sound, clip_idx) = sdr_rec - sdr_clip;


            disp(['Done: ' num2str(cnt), '/', num2str(numel(SDR))]);
        end
        save('results/speech_test_references/CSL1_Sounds.mat')

        dSDR_all_csl1 = dSDR;
        save('survey toolbox/Numerical_results/dSDR_all_declippingResults_speech.mat', 'dSDR_all_csl1', '-append')
    end
end