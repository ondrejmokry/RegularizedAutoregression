%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        DEQUANTIZATION          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        WHOLE DATABASE          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pavel Záviška, Brno University of Technology, 2020
%
% dequantization_whole_database.m adapted for the use of the regularized AR
% model by Ondøej Mokrý, Brno University of Technology, 2025
%
% To reproduce the dSDR values from paper, enable the dsdr_decterm 
% parameter, which terminates the algorithms after a SDR drop. On the other
% hand, to reproduce the PEMO-Q ODG values, let this parameter disabled and
% run full 500 iterations. [https://github.com/zawi01/audio_dequantization]

clear
clc
close all
%#ok<*UNRCH>

% using toolbox LTFAT
% ltfatstart

filename = 'dequant_all';
directory = ['results/', filename];
[~, ~] = mkdir(directory);

addpath('dequantization toolbox/Algorithms')
addpath('dequantization toolbox/Sounds')
addpath('dequantization toolbox/Tools')
addpath('utils')
addpath(genpath('survey toolbox/Evaluation_algorithms/PemoQ'))

algorithms = { ...
    'DR_cons_l1_syn', ...
    'CP_cons_l1_ana', ...
    'A_SPADQ', ...
    'S_SPADQ', ...
    'S_SPADQ_DR', ...
    'FISTA_incons_l1_syn', ...
    'DR_incons_l1_syn', ...
    'CP_incons_l1_ana', ...
    'DR_incons_l1_ana', ...
    'FISTA_incons_l1_ana', ...
    'Janssen_cons', ...
    'Janssen_incons', ...
    'Janssen_sparse_cons', ...
    'Janssen_sparse_incons'};

sounds = { ...
    'a08_violin', ...
    'a16_clarinet', ...
    'a18_bassoon', ...
    'a25_harp', ...
    'a35_glockenspiel', ...
    'a41_celesta', ...
    'a42_accordion', ...
    'a58_guitar_sarasate', ...
    'a60_piano_schubert', ...
    'a66_wind_ensemble_stravinsky'};

% here set algorithms, sounds, and wordlenghts to compute
alg_idxs = 1:14;
sound_idxs = 1:length(sounds);      
wordlengths = 3:7;    

STORE_dSDR_PROCESS = true;
STORE_OBJ_PROCESS = true;
STORE_DEQ_SOUNDS = true;

% initialization of counter
cnt = 0;
cases = length(alg_idxs)*length(sound_idxs)*length(wordlengths);

for algorithm = alg_idxs

    % initialization of matrices
    SDR_mat = NaN(length(sounds), length(wordlengths));
    ODG_mat = NaN(length(sounds), length(wordlengths));
    dSDR_mat = NaN(length(sounds), length(wordlengths));
    dODG_mat = NaN(length(sounds), length(wordlengths));
    time_mat = NaN(length(sounds), length(wordlengths));

    if STORE_dSDR_PROCESS
        dSDR_process_mat = cell(length(sounds), length(wordlengths));
    end
    
    if STORE_OBJ_PROCESS
        objective_process_mat = cell(length(sounds), length(wordlengths));
    end

    if isfile(['results/', filename, '.mat'])
        load(['results/', filename, '.mat'])
        if isfield(SDR, algorithms{algorithm})
            SDR_mat = SDR.(algorithms{algorithm});
            dSDR_mat = dSDR.(algorithms{algorithm});
            ODG_mat = ODG.(algorithms{algorithm});
            dODG_mat = dODG.(algorithms{algorithm});
            time_mat = TIME.(algorithms{algorithm});
        end
    end

    for sound = sound_idxs
        for wl = 1:length(wordlengths)
            reconstructions = struct;

            cnt = cnt + 1;
            if ~isnan(SDR_mat(sound, wl))
                continue
            end

            %% input file settings
            [data, fs] = audioread("dequantization toolbox\Sounds\" + sounds{sound} + ".wav");

            % peak-normalization
            maxAbsVal = max(abs(data));
            data = data/maxAbsVal;

            % signal length
            param.Ls = length(data);


            %% General settings
            param.wordlength = wordlengths(wl);     % set the wordlength in bits
            param.algorithm = algorithms{algorithm}; % algorithm to compute declipping, options: 'DR', 'CP', 'A-SPADQ', 'S-SPADQ', 'S-SPADQ_DR', 'FISTA'


            %% Settings for l1-minimization algorithms (CP, DR)
            if any(strcmp(param.algorithm, {'DR_cons_l1_syn', 'CP_cons_l1_ana', 'FISTA_incons_l1_syn', ...
                'DR_incons_l1_syn', 'CP_incons_l1_ana', 'DR_incons_l1_ana', 'FISTA_incons_l1_ana'}))
               
                % frame settings
                param.w = 8192;
                param.a = param.w/4;
                param.M = 2*param.w; % M >= w
                param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

                % construction of frame
                param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M));
                param.F = frameaccel(param.F, param.Ls);  % precomputation for a fixed signal length

                % general settings of the l1 minimization algorithms (algorithm parameters are set directly in the respective m-file.)
                paramsolver.maxit = 500;    % maximum number of iterations
                paramsolver.minit = 25;    % minimum number of iterations 
                paramsolver.verbose = 0;    % display parameter
                paramsolver.comp_dsdr = STORE_dSDR_PROCESS;  % compute and store dSDR during iterations
                paramsolver.dsdr_decterm = 1;  % terminate algorithm if the SDR value starts to decrease
                paramsolver.comp_obj = STORE_OBJ_PROCESS;   % compute and store objective function values during iterations
    
            end

            %% Settings for SPADQ algorithms
            if any(strcmp(param.algorithm, {'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR'}))
            
                % window parameters
                param.w = 8192;       % window length
                param.a = param.w/4;  % window shift
                param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

                % DFT parameters
                param.F = frame('dft');
                param.F.redundancy = 2;  %non-native, our own additional parameter
                param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
                param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);
                
                % general settings of the SPADQ algorithms
                paramsolver.verbose = 0;
                paramsolver.comp_sdr = STORE_dSDR_PROCESS;
                paramsolver.comp_obj = STORE_OBJ_PROCESS;
            end

            %% Settings for Janssen
            if any(strcmp(param.algorithm, {'Janssen_cons', 'Janssen_incons', ...
                    'Janssen_sparse_cons', 'Janssen_sparse_incons'}))

                % window parameters
                param.w = 2048;       % window length
                param.a = param.w/4;  % window shift
                param.wtype = 'rect'; % probably better for Janssen
                param.segmentation = true;
 
                % AR model parameters
                param.p = 1024;

                % algorithm parameters
                param.iterations = 2;
                param.coefextra = false;
                param.sigextra = false;
                param.coefaccel = true;
                param.sigaccel = true;
                param.linesearch = false;
                param.plotLS = false;
                param.gammaC = 0.1;
                param.gammaS = 10;
                param.DRmaxit = 1000;
                param.mat = "toeplitz";

                % regularization parameters
                param.lambdaC = 0.1;
                param.lambdaS = 1.0;

            end

            %% Quantization
            [data_quant, param.delta] = quant(data, param.wordlength); % quantizing the original signal and computing the quantization step


            %% Optimization algorithm

            timer = tic;

            switch param.algorithm
                % consistent l1-minimization using synthesis model of the signal, Douglas-Rachford algorithm      
                case {'DR_cons_l1_syn'}
                    [data_rec, dsdr_iter, obj_iter] = dr_cons_l1_syn(data_quant, param, paramsolver, data);

                % consistent l1-minimization using analysis model of the signal, Chambolle-Pock algorithm
                case {'CP_cons_l1_ana'}
                    [data_rec, dsdr_iter, obj_iter] = cp_cons_l1_ana(data_quant, param, paramsolver, data);

                % non-convex l0-minimization based on ADMM, SPADQ algorithms
                case {'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR'}
                    % paramsolver parameters
                    paramsolver.s = 1;   % increment of k
                    paramsolver.r = 1;   % every r-th iteration increment k by s
                    paramsolver.epsilon = 0.01;  % stopping criterion of termination function
                    paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations

                    [data_rec, dsdr_iter, obj_iter] = spadq_segmentation(data_quant, param, paramsolver, data);

                % inconsistent l1-minimization using synthesis model of the signal, FISTA
                case {'FISTA_incons_l1_syn'}
                    [data_rec, dsdr_iter, obj_iter] = fista_incons_l1_syn(data_quant, param, paramsolver, data);

                % inconsistent l1-minimization using synthesis model of the signal, Douglas-Rachford algorithm
                case {'DR_incons_l1_syn'}
                    [data_rec, dsdr_iter, obj_iter] = dr_incons_l1_syn(data_quant, param, paramsolver, data);

                % inconsistent l1-minimization using analysis model of the signal, Chambolle-Pock algorithm   
                case {'CP_incons_l1_ana'}
                    [data_rec, dsdr_iter, obj_iter] = cp_incons_l1_ana(data_quant, param, paramsolver, data);        

                % inconsistent l1-minimization using analysis model of the signal, Douglas-Rachford algorithm
                case {'DR_incons_l1_ana'}
                    [data_rec, dsdr_iter, obj_iter] = dr_incons_l1_ana(data_quant, param, paramsolver, data);

                % inconsistent l1-minimization using analysis model of the signal, FISTA
                case {'FISTA_incons_l1_ana'}
                    [data_rec, dsdr_iter, obj_iter] = fista_incons_l1_ana(data_quant, param, paramsolver, data);

                % regularized autoregression
                case {'Janssen_cons', 'Janssen_incons', 'Janssen_sparse_cons', 'Janssen_sparse_incons'}
                    masks = struct;
                    if contains(param.algorithm, 'incons')
                        lambdaS = param.lambdaS;
                    else
                        lambdaS = Inf;
                    end
                    if contains(param.algorithm, 'sparse')
                        lambdaC = param.lambdaC;
                    else
                        lambdaC = 0;
                    end

                    data_rec = janssen("dequantization", ...
                        data_quant, masks, [lambdaC, lambdaS], param.p, param.iterations, ...
                        "delta", param.delta, "segmentation", param.segmentation, ...
                        "a", param.a, "w", param.w, "wtype", param.wtype, ...
                        "DRmaxit", param.DRmaxit, "decompose", true, "mat", param.mat, ...
                        "gammaC", param.gammaC, "gammaS", param.gammaS, ...
                        "coefaccel", param.coefaccel, "sigaccel", param.sigaccel, ...
                        "coefextra", param.coefextra, "sigextra", param.sigextra, ...
                        "linesearch", param.linesearch, "plotLS", param.plotLS, ...
                        "saveall", true, "verbose", false);

                    % this is computationally too expensive
                    obj_iter = [];

                otherwise
                    error('Invalid algorithm is set!');
            end

            time = toc(timer);
            
            %% Time & SDR evaluation

            % compute dSDR in iterations for AR methods
            if contains(param.algorithm, 'Janssen')
                dsdr_iter = zeros(size(data_rec, 2), 1);
                sdr_quant = snr(data, data_quant-data);
                for i = 1:size(data_rec, 2)
                    dsdr_iter(i) = snr(data, data_rec(:, i)-data) - sdr_quant;
                end
                data_rec = data_rec(:, end);
            end
            
            % rename and save dequantized sound
            if STORE_DEQ_SOUNDS
                reconstructions.([sounds{sound} '_0' num2str(wordlengths(wl)) '_rec_' algorithms{algorithm}]) = data_rec;
            end
            
            % store dSDR course through iterations
            if STORE_dSDR_PROCESS, dSDR_process_mat{sound, wl} = dsdr_iter; end
            
            % store course of the objective function through iterations
            if STORE_OBJ_PROCESS, objective_process_mat{sound, wl} = obj_iter; end
            
            % store computational time
            time_mat(sound, wl) = time;
            
            % compute and store the SDR and dSDR values of the
            % reconstructed (dequantized) signal
            sdr_quant = sdr(data, data_quant);
            sdr_rec = sdr(data, data_rec);
            SDR_mat(sound, wl) = sdr_rec;
            dSDR_mat(sound, wl) = sdr_rec - sdr_quant;

            % compute ODG
            [~, ~, odg_quant, ~] = audioqual_silent(data, data_quant, fs);
            [~, ~, odg_rec, ~] = audioqual_silent(data, data_rec, fs);
            ODG_mat(sound, wl) = odg_rec;
            dODG_mat(sound, wl) = odg_rec - odg_quant;
            
            disp(['Done: ', num2str(cnt), ' / ', num2str(cases)])

            SDR.(algorithms{algorithm}) = SDR_mat;
            dSDR.(algorithms{algorithm}) = dSDR_mat;
            ODG.(algorithms{algorithm}) = ODG_mat;
            dODG.(algorithms{algorithm}) = dODG_mat;
            TIME.(algorithms{algorithm}) = time_mat;

            save(['results/', filename, '.mat'], ...
                'algorithms', 'sounds', 'wordlengths', ...
                'SDR', 'dSDR', 'ODG', 'dODG', 'TIME')
            
            if STORE_dSDR_PROCESS
                dSDR_process.(algorithms{algorithm}) = dSDR_process_mat;
                save(['results/', filename, '.mat'], 'dSDR_process', '-append')
            end
            
            if STORE_OBJ_PROCESS
                objective_process.(algorithms{algorithm}) = objective_process_mat;
                save(['results/', filename, '.mat'], 'objective_process', '-append')
            end

            if STORE_DEQ_SOUNDS
                if isfile([directory '/' sounds{sound}(1:3) '_0' num2str(wordlengths(wl)) '.mat'])
                    save([directory '/' sounds{sound}(1:3) '_0' num2str(wordlengths(wl)) '.mat'], '-struct', 'reconstructions', '-append')
                else
                    save([directory '/' sounds{sound}(1:3) '_0' num2str(wordlengths(wl)) '.mat'], '-struct', 'reconstructions')
                end
            end

        end
    end
 
end