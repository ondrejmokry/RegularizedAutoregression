%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%         Social Sparsity        %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%           (all sounds)         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation of the Social Sparsity Declipper was kindly provided
% by Matthieu Kowalski.
% The codes were then slightly edited to fit the common interface of the
% declipping toolbox by
% Pavel Záviška, Brno University of Technology, 2020
% adapted by Ondøej Mokrý for the speech dataset, 2025
% 
% using toolbox LTFAT
ltfatstart;

addpath('survey toolbox/Tools')
addpath(genpath('survey toolbox/Evaluation_algorithms'))
addpath('survey toolbox/Methods/Social Sparsity')

close all
clear variables

% Database
sounds = dir('speech');
sounds = sounds(3:end-1);

shrinkage = 'PEW'; % 'L', 'WGL', 'EW', 'PEW';
input_SDRs = [ 3, 5, 7, 10 ];

verbose = false;

% initialization of matrices for SDR and dSDR values, and computational time
SDR = NaN(length(sounds), length(input_SDRs));
dSDR = NaN(length(sounds), length(input_SDRs));
TIME = NaN(length(sounds), length(input_SDRs));

% initialization of counter
cnt = 0;


for sound = 1:length(sounds)
    for clip_idx = 1:length(input_SDRs)
        
        cnt = cnt + 1;
        % load audio-file
        [data, fs] = audioread(['speech/' sounds(sound).name]);
        
        % signal length
        Ls = length(data);
        
        
        %% settings
        % input SDR of the clipped signal
        inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold
        
        % DGT parameters
        wtype = 'hann';   % window type
        w = 2048;         % window length
        a = w / 4;  % window shift
        M = 2*2048;       % number of frequency channels
        
        % clipping
        [data_clipped, masks, theta, ~, ~] = clip_sdr(data, inputSDR); % clipping of original signal
        
        data_r = zeros(Ls,1);
        data_r(masks.Mr) = data_clipped(masks.Mr);
        
        data_c = zeros(Ls,1);
        data_c(masks.Mh) = theta;
        data_c(masks.Ml) = -theta;
        
        % construction of the frame
        g = gabwin({'tight', {'hann', w}}, a, M);
        
        G_clip = dgtreal(data_clipped, g, a, M);
        [nf,nt] = size(G_clip);
        
        gab.analysis = @(x) dgtreal(x,g,a,M);
        gab.synthesis = @(x) idgtreal(x,g,a,M,Ls);
        
        % Setting of the Shrinkage operator
        shrinkL   = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', 1, 'center', [1 1], 'expo', 1,'orth',0);
        shrinkWGL = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', ones(3,7), 'center', [2,4], 'expo', 1,'orth',0);
        shrinkEW  = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', 1, 'center', [1 1], 'expo', 2,'orth',0);
        shrinkPEW = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', ones(3,7), 'center', [2,4], 'expo', 2,'orth',0);
        
        switch shrinkage
            case 'L'
                shrink = shrinkL;
            case 'WGL'
                shrink = shrinkWGL;
            case 'EW'
                shrink = shrinkEW;
            case 'PEW'
                shrink = shrinkPEW;
            otherwise
                error('Invalid shrinkage operator!')
        end
        
        Gdeclip = zeros(nf,nt);
        GdeclipZ = Gdeclip;
        
        %% ISTA
        tic
        
        number_lambdas = 20;
        inner_iterations = 500;
        delta = 1e-4;
        iter_cnt = 0;
        
        dSDR_process = NaN(number_lambdas*inner_iterations, 1);
        data_rec = data_clipped;
        
        for lambda=logspace(-1,-4,number_lambdas)
            shrink.lambda = lambda;
            for k=1:inner_iterations
                
                data_rec_old = data_rec;
                iter_cnt = iter_cnt + 1;
                
                % forward step
                GdeclipOLD = Gdeclip;
                
                xdeclipZ = gab.synthesis(GdeclipZ);
                
                r1 = zeros(Ls,1);
                r1(masks.Mr) = data_r(masks.Mr) - xdeclipZ(masks.Mr);
                
                grad1 = -gab.analysis(r1);
                
                r2 = zeros(Ls,1);
                r2(masks.Mh) = data_c(masks.Mh) - xdeclipZ(masks.Mh);
                r2(masks.Ml) = data_c(masks.Ml) - xdeclipZ(masks.Ml);
                
                r2(abs(xdeclipZ)>theta) = 0;
                grad2 = -gab.analysis(r2);
                
                
                GdeclipZ = GdeclipZ - grad1 - grad2;
                
                
                % thresholding step
                Gdeclip = gen_thresh(GdeclipZ, shrink);
                
                
                % relaxation step
                GdeclipZ = Gdeclip + (k-1)/(k+5) * (Gdeclip - GdeclipOLD);
                
                
                
                data_rec = gab.synthesis(Gdeclip);
                if norm(data_rec_old - data_rec) < delta
                    break
                end
                
                                
                if verbose
                    sdr_rec = sdr(data, data_rec); %#ok<UNRCH>
                    fprintf('  lambda = %f -- it = %d -- SDR = %f -- dSDR = %f \n',lambda,k,sdr_rec,sdr_rec-input_SDRs(clip_idx));
                end
                
            end
        end
        
        time = toc;
        
        %% Evaluation of the results
        % Time
        TIME(sound, clip_idx) = time;
        
        % Save restored file
        if input_SDRs(clip_idx) < 10
            eval([sounds(sound).name(1:end-4) '_rec_SS_' shrinkage '_0' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
        else
            eval([sounds(sound).name(1:end-4) '_rec_SS_' shrinkage '_' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
        end
        
        % SDR
        sdr_clip = sdr(data, data_clipped);
        sdr_rec = sdr(data, data_rec);
        
        SDR(sound, clip_idx) = sdr_rec;
        dSDR(sound, clip_idx) = sdr_rec - sdr_clip;
        
        
        disp(['Done: ' num2str(cnt), '/', num2str(numel(SDR))]);
        
        save('results/speech_test_references/SS_PEW_Sounds_ALL.mat');

        dSDR_all_SS_PEW = dSDR;
        save('survey toolbox/Numerical_results/dSDR_all_declippingResults_speech.mat', 'dSDR_all_SS_PEW', '-append')
        
    end
end