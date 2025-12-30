clear
clc
close all

% name of the file to load/save
filename = 'speech_test';

% name of the directory to load the waveforms from
directory = ['results/', filename];

%% load matrics and hyperparameters
R = load(['results/', filename, '.mat']);

%% prepare the new fields to R.reconstructed
for algo = 1:length(R.algos)
    if algo == 2
        dims = [length(R.audio_files), length(R.input_SDRs), length(R.lambdaC), length(R.lambdaS)];
    else
        dims = [length(R.audio_files), length(R.input_SDRs), length(R.lambdaC)];
    end
    R.reconstructed(algo).STOIs = NaN(dims);
    R.reconstructed(algo).MOSs = NaN(dims);
    R.reconstructed(algo).NSIMs = NaN(dims);
    R.reconstructed(algo).SDRs_check = NaN(dims);
end

%% iterate over files and solutions
for i = 1:length(R.audio_files)
    [data, fs] = audioread(['speech/', R.audio_files(i).name]);
    for j = 1:length(R.input_SDRs)
        S = load([directory, '/', R.audio_files(i).name(1:end-4), '_', num2str(R.input_SDRs(j), '%02d')]);
        for m = 1:length(R.lambdaC)
            for n = 1:length(R.lambdaS)
                for algo = 1:length(R.algos)
            
                    % skip lambdaS iteration for inpainting and glp
                    if algo ~= 2 && n > 1
                        continue
                    end

                    % parse the file
                    switch algo
                        case 1, declipped = S.(['inp_', num2str(10*m)]);
                        case 2, declipped = S.(['dec_', num2str(10*m + n)]);
                        case 3, declipped = S.(['glp_', num2str(10*m)]);
                    end

                    % evaluate
                    R.reconstructed(algo).STOIs(i, j, m, n) = stoi(declipped, data, fs);
                    mn = visqol(declipped, data, fs, 'Mode', 'speech', 'OutputMetric', 'MOS and NSIM');
                    R.reconstructed(algo).MOSs(i, j, m, n) = mn(1);
                    R.reconstructed(algo).NSIMs(i, j, m, n) = mn(2);
                    R.reconstructed(algo).SDRs_check(i, j, m, n) = snr(data, data-declipped);
                end
            end
        end
    end
end
save(['results/', filename, '.mat'], '-struct', 'R')

%% reference algos
refalgos = {'ASPADE', 'CSL1', 'SS_PEW'};
refvars = {'aspade', 'csl1', 'SS_PEW'};

STOI = struct;
MOS = struct;
NSIM = struct;
dSDR_check = struct;

for i = 1:length(refalgos)
    
    % load audio from here
    refdata = load(['results/speech_test_references/', refalgos{i}, '_Sounds.mat']);

    % iterate over files
    STOI.(['STOI_', refvars{i}]) = zeros(length(R.audio_files), length(R.input_SDRs));
    MOS.(['MOS_', refvars{i}]) = zeros(length(R.audio_files), length(R.input_SDRs));
    NSIM.(['NSIM_', refvars{i}]) = zeros(length(R.audio_files), length(R.input_SDRs));
    dSDR_check.(['dSDR_all_', refvars{i}]) = zeros(length(R.audio_files), length(R.input_SDRs));

    for j = 1:length(R.audio_files)
        [data, fs] = audioread(['speech/', R.audio_files(j).name]);
        for k = 1:length(R.input_SDRs)
            % load declipped
            declipped = refdata.([R.audio_files(j).name(1:end-4), '_rec_', refvars{i}, '_', num2str(R.input_SDRs(k), '%02d')]);

            % evaluate
            STOI.(['STOI_', refvars{i}])(j, k) = stoi(declipped, data, fs);
            mn = visqol(declipped, data, fs, 'Mode', 'speech', 'OutputMetric', 'MOS and NSIM');
            MOS.(['MOS_', refvars{i}])(j, k) = mn(1);
            NSIM.(['NSIM_', refvars{i}])(j, k) = mn(2);
            % dSDR_check.(['dSDR_all_', refvars{i}])(j, k) = snr(data, declipped-data) - R.input_SDRs(k);
        end
    end

    % save
    save('survey toolbox/Numerical_results/STOI_declippingResults_speech.mat', ...
        '-struct', 'STOI');
    save('survey toolbox/Numerical_results/MOS_declippingResults_speech.mat', ...
        '-struct', 'MOS');
    save('survey toolbox/Numerical_results/NSIM_declippingResults_speech.mat', ...
        '-struct', 'NSIM');
    % save('survey toolbox/Numerical_results/dSDR_check_declippingResults_speech.mat', ...
    %     '-struct', 'dSDR_check');

end