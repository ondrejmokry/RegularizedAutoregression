%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%              SPADE             %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%          (all sounds)          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pavel Záviška, Brno University of Technology, 2020
% adapted by Ondøej Mokrý for the speech dataset, 2025
% 
% using toolbox LTFAT
ltfatstart

addpath('survey toolbox/Tools')
addpath(genpath('survey toolbox/Evaluation_algorithms'))
addpath('survey toolbox/Methods/SPADE')

close all
clear variables

% Database
sounds = dir('speech');
sounds = sounds(3:end-1);

algorithms = {'aspade', 'sspade_new'};
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
            for redundancy = 2
                
                cnt = cnt+1;
                %% input file
                % load audio-file              
                [data, fs] = audioread(['speech/' sounds(sound).name]);

                % signal length
                param.Ls = length(data);
                
                
                %% settings
                param.inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold
                
                % window parameters
                param.w = 2048 ;       % window length
                param.a = param.w/4;  % window shift
                param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html
                
                % DFT parameters
                param.F = frame('dft');
                param.F.redundancy = redundancy;  %non-native, our own additional parameter
                param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
                param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);
                               
                % paramsolver parameters
                paramsolver.s = 1;   % increment of k
                paramsolver.r = 2;   % every r-th iteration increment k by s   
                paramsolver.epsilon = 0.1;  % stopping criterion of termination function
                paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations

                paramsolver.store_sdr = 0; 
                paramsolver.store_obj = 0;
                
                paramsolver.verbose = false; % false - no listing, true - inform about processed blocks
                
                % algorithm
                param.algorithm = algorithms{algorithm}; % algorithm to compute declipping, options: 'aspade', 'sspade', 'sspade_new'
                
 
                %% clipping
                [data_clipped, param.masks, param.theta, ~, ~] = clip_sdr(data, param.inputSDR); % clipping of original signal

                
                %% Main algorithm            
                tic;
                
                [data_rec, ~, ~] = spade_segmentation(data_clipped, param, paramsolver, data);
                
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

            if strcmp(algorithms{algorithm}, 'aspade')
                save('results/speech_test_references/ASPADE_Sounds.mat');
                dSDR_all_aspade = dSDR;
                save('survey toolbox/Numerical_results/dSDR_all_declippingResults_speech.mat', 'dSDR_all_aspade', '-append')
            else
                save('results/speech_test_references/SSPADE_DR_Sounds.mat');
                dSDR_all_sspade_new = dSDR;
                save('survey toolbox/Numerical_results/dSDR_all_declippingResults_speech.mat', 'dSDR_all_sspade_new', '-append')
            end
        end        
    end
end