clear
clc
close all

fold = 'results/survey_test_rect';
toolboxpath = 'survey toolbox';
addpath(genpath(toolboxpath))
addpath('utils')

% load metrics
metrics = load([fold, '.mat']);

% setup
audio_files = { 'a08_violin', ...
                'a16_clarinet', ...
                'a18_bassoon', ...
                'a25_harp', ...
                'a35_glockenspiel', ...
                'a41_celesta', ...
                'a42_accordion', ...
                'a58_guitar_sarasate', ...
                'a60_piano_schubert', ...
                'a66_wind_ensemble_stravinsky' };
nofiles = length(audio_files);
input_SDRs = [ 5, 7, 10, 15 ];
tol = 1e-8;

% prepare algo names and shortcuts
algos = cell(1, 12); 
algonames = cell(1, 12);
counter = 0;
for i = 1:length(metrics.lambdaS)
    for j = 1:length(metrics.lambdaC)
        counter = counter + 1;
        algos{counter} = ['dec_', num2str(10*j + i)];
        algonames{counter} = ['dec., $\lambda_C$ = ', num2str(metrics.lambdaC(j)), ', $\lambda_S$ = ', num2str(metrics.lambdaS(i))];
    end
end
for i = 1:length(metrics.lambdaC)
    counter = counter + 1;
    algos{counter} = ['glp_', num2str(10*i)];
    algonames{counter} = ['GLP, $\lambda_C$ = ', num2str(metrics.lambdaC(i))];
end
for i = 1:length(metrics.lambdaC)
    counter = counter + 1;
    algos{counter} = ['inp_', num2str(10*i)];
    algonames{counter} = ['inp., $\lambda_C$ = ', num2str(metrics.lambdaC(i))];
end

% reorder
algos = algos([10:12, 7:9, 1:6]);
algonames = algonames([10:12, 7:9, 1:6]);

% prepare tables
R = table('Size', [length(audio_files)*length(input_SDRs), 2+length(algos)], ...
    'VariableTypes', ['string', 'int8', repmat({'double'}, 1, length(algos))], ...
    'VariableNames', ['audio_file', 'input_SDR', algos]);
H = R;
L = R;
N = R;
Nlog = R;
SDR = R;
PEMOQ = R;
PEAQ = R;

row = 1;
for i = 1:length(audio_files)

    %% load audio file
    fprintf(['Processing audio ''', audio_files{i} , '.wav''\n']);
    [data, fs] = audioread([toolboxpath, '/Sounds/' audio_files{i} '.wav']);

    for j = 1:length(input_SDRs)

        R.audio_file{row} = audio_files{i};
        R.input_SDR(row) = input_SDRs(j);
        H.audio_file{row} = audio_files{i};
        H.input_SDR(row) = input_SDRs(j);
        L.audio_file{row} = audio_files{i};
        L.input_SDR(row) = input_SDRs(j);
        N.audio_file{row} = audio_files{i};
        N.input_SDR(row) = input_SDRs(j);
        Nlog.audio_file{row} = audio_files{i};
        Nlog.input_SDR(row) = input_SDRs(j);
        SDR.audio_file{row} = audio_files{i};
        SDR.input_SDR(row) = input_SDRs(j);
        PEMOQ.audio_file{row} = audio_files{i};
        PEMOQ.input_SDR(row) = input_SDRs(j);
        PEAQ.audio_file{row} = audio_files{i};
        PEAQ.input_SDR(row) = input_SDRs(j);

        %% clip audio file (using the function in survey toolbox/Tools)
        [data_clipped, masks, theta, trueSDR, percentage] = clip_sdr(data, input_SDRs(j));
        
        %% load reconstruction
        declipped = load(sprintf('%s/%s_%02d.mat', fold, audio_files{i}(1:3), input_SDRs(j)));

        %% evaluate consistency
        for k = 1:length(algos)

            algo = algos{k};

            reco = declipped.(algo);
            
            R.(algo)(row) = sum(abs(reco(masks.Mr) - data(masks.Mr)) < tol)/sum(masks.Mr);
            H.(algo)(row) = sum(reco(masks.Mh) >= theta-tol)/sum(masks.Mh);
            L.(algo)(row) = sum(reco(masks.Ml) <= -theta+tol)/sum(masks.Ml);
            
            SDR.(algo)(row) = snr(data, data-reco);

            proj = reco;
            proj(masks.Mr) = data(masks.Mr);
            proj(masks.Mh) = max(theta, proj(masks.Mh));
            proj(masks.Ml) = min(-theta, proj(masks.Ml));
            N.(algo)(row) = max(0.5*norm(proj - reco)^2, 1e-12);
            Nlog.(algo)(row) = log10(0.5*norm(proj - reco)^2);

        end
        row = row + 1;
        
    end
end

%% plot SDR versus consistency

% perm = [7:12, 4:6, 1:3];
% load('utils/barcolors.mat')
% colors = colors(perm, :);

colors = turbo(12);

figure
scatter(SDR{:,3:end} - double(SDR.input_SDR(:)), N{:, 3:end}, [], colors, 'filled')
xlabel('dSDR')
ylabel('projection distance')
set(gca, 'YScale', 'log')
grid on
legend(algonames)
box on

figure
b = boxchart(Nlog, 3:14);
for i = 1:length(b)
    b(i).BoxFaceColor = colors(i, :);
end
set(gca, 'XTickLabel', algonames)
box on
grid on
ylabel('log projection distance')
legend(algonames, Location='southeast')