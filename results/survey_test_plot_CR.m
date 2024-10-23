% the script plots the data generated by survey_test.m, 
% such that both "CR" and "no CR" options are plotted using stacked bar plot 

clear
clc
close all
%#ok<*UNRCH>
addpath('../utils')

plotmetrics = 1:3; % 1: SDR, 2: PEMO-Q ODG, 3: PEAQ ODG

% the audio signals to choose from
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
           
%% load the data
filename = 'survey_test_rect';
A = load(filename);
B = load([filename, '_CR']);
almostzero = 1e-4;

%% reference algorithms
% available options
dec_algos = { 'consOMP', ...       % 1
              'aspade', ...        % 2
              'sspade_new', ...    % 3
              'CP', ...            % 4
              'DR', ...            % 5
              'SS_EW', ...         % 6
              'SS_PEW', ...        % 7
              'csl1', ...          % 8
              'pcsl1', ...         % 9
              'pwcsl1', ...        % 10
              'reweighted_CP', ... % 11
              'reweighted_DR', ... % 12
              'CP_parabola', ...   % 13
              'DR_parabola', ...   % 14
              'DL', ...            % 15
              'NMF'...            % 16
              };

% per each input SDR, choose the reference algorithms
% references = { [ 2, 3, 7 ], ... % 5 dB
%                [ 2, 7, 16 ], ... % 7 dB
%                [ 7, 13, 16 ], ... % 10 dB
%                [ 7, 13, 16 ], ... % 15 dB
%                };

references = 1:16;

% find the unique references
if iscell(references)
    unique_refs = unique([references{:}]);
else
    unique_refs = unique(references(:));
end

% save some useful numbers
num_inp = length(A.lambdaC);
num_glp = length(A.lambdaC);
num_dec = length(A.lambdaC)*length(A.lambdaS);
num_ref = length(unique_refs);
num_tot = num_inp + num_glp + num_dec + num_ref;

%% prepare the legend entries
legstr = cell(num_tot, 1);
counter = 0;

% inpainting
for i = 1:num_inp
    counter = counter + 1;
    legstr{counter} = ['inpainting, $$\lambda_\textup{C}$$ = ', num2str(A.lambdaC(i))];
end

% glp
for i = 1:num_glp
    counter = counter + 1;
    legstr{counter} = ['GLP, $$\lambda_\textup{C}$$ = ', num2str(A.lambdaC(i))];
end

% declipping
for i = 1:length(A.lambdaS)
    for j = 1:length(A.lambdaC)
        counter = counter + 1;
        legstr{counter} = ['declipping, $$\lambda_\textup{C}$$ = ', num2str(A.lambdaC(j)), ', $$\lambda_\textup{S}$$ = ', num2str(A.lambdaS(i))];
    end
end

% references
for i = 1:num_ref
    counter = counter + 1;
    legstr{counter} = strrep(dec_algos{unique_refs(i)}, '_', ' ');
end

%% prepare colors
% set colormap
% mapcolors = turbo(1 + length(A.lambdaS) + 1 + num_ref);
mapcolors = turbo(num_ref);
mapcolors = [mapcolors(1:length(A.lambdaS)+2, :); mapcolors];
colors = zeros(num_tot, 3);
colors(num_inp + num_glp + num_dec + 1:end, :) = mapcolors(length(A.lambdaS)+3:end, :);
for i = 1:length(A.lambdaC)
    frac = i/length(A.lambdaC);

    % inpainting
    colors(i, :) = mapcolors(1, :)*frac + [1 1 1]*(1-frac);

    % glp
    colors(num_inp + i, :) = mapcolors(2, :)*frac + [1 1 1]*(1-frac);

    % declipping
    for j = 1:length(A.lambdaS)
        colors(num_inp + num_glp + (j-1)*length(A.lambdaC) + i, :) = mapcolors(2+j, :)*frac + [1 1 1]*(1-frac);
    end
end

%% plot mean from all signals
fcnt = 0; % figure counter
for metric = plotmetrics

    % initialize plotdata (the mean values)
    plotdata = NaN(length(A.input_SDRs), num_tot, 2);

    % switch metric (what a surprise)
    switch metric
        case 1, field = 'SDRs';   reflabel = 'dSDR_all'; base = 0;
        case 2, field = 'PEMOQs'; reflabel = 'PEMO_Q';   base = -4;
        case 3, field = 'PEAQs';  reflabel = 'PEAQ';     base = -4;
    end

    % read the inpainting + glp + declipping data
    for method = 1:3

        % read the data
        switch method
            case 1
                data = A.reconstructed(1).(field);
                dataCR = B.reconstructed(1).(field);
            case 2
                data = A.reconstructed(3).(field);
                dataCR = B.reconstructed(3).(field);
            case 3
                data = A.reconstructed(2).(field);
                dataCR = B.reconstructed(2).(field);
        end

        % compute dSDR (if applicable)
        if metric == 1
            if method == 3
                data = data - repmat(A.clipped.SDRs, 1, 1, length(A.lambdaC), length(A.lambdaS));
                dataCR = dataCR - repmat(A.clipped.SDRs, 1, 1, length(A.lambdaC), length(A.lambdaS));
            else
                data = data - repmat(A.clipped.SDRs, 1, 1, length(A.lambdaC));
                dataCR = dataCR - repmat(A.clipped.SDRs, 1, 1, length(A.lambdaC));
            end
        end
        data = mean(data, 1);
        data = squeeze(data);
        if metric > 1
            dataCR = mean(dataCR, 1);
            dataCR = squeeze(dataCR) + almostzero; % it causes problems if the two are equal
        else
            dataCR = data;
        end

        % save to plotdata
        switch method
            case 1
                for i = 1:num_inp
                    plotdata(:, i, 1) = data(:, i);
                    plotdata(:, i, 2) = dataCR(:, i);
                end
            case 2
                for i = 1:num_glp
                    plotdata(:, num_inp + i, 1) = data(:, i);
                    plotdata(:, num_inp + i, 2) = dataCR(:, i);
                end
            case 3
                for i = 1:length(A.lambdaS)
                    for j = 1:length(A.lambdaC)
                        plotdata(:, num_inp + num_glp + (i-1)*length(A.lambdaC) + j, 1)...
                            = data(:, j, i);
                        plotdata(:, num_inp + num_glp + (i-1)*length(A.lambdaC) + j, 2)...
                            = dataCR(:, j, i);  
                    end
                end
        end
    end

    % read the reference data
    for i = 1:length(A.input_SDRs)
        for j = 1:num_ref
            refdata = load(['../survey toolbox/Numerical_results/', reflabel, '_declippingResults.mat']);
            refdata = refdata.([reflabel, '_', dec_algos{unique_refs(j)}]);
            refdata = mean(refdata, 1);
            
            if metric > 1
                % these data are available only for PEAQ, PEMO-Q and SDR on
                % clipped samples, but we compute SDR on the whole signal
                refdataCR = load(['../survey toolbox/Numerical_results/results_RR_web/', reflabel, '_declippingResults_CR.mat']);
                refdataCR = refdataCR.([reflabel, '_', dec_algos{unique_refs(j)}, '_CR']);
                refdataCR = mean(refdataCR, 1) + almostzero; % it causes problems if the two are equal
            else
                refdataCR = refdata;
            end
            
            % the index for refdata needs to be shifted, because the input
            % SDRs of 1 dB and 3 dB were omitted in the test
            plotdata(i, num_inp + num_glp + num_dec + j, 1) = refdata(2+i);
            plotdata(i, num_inp + num_glp + num_dec + j, 2) = refdataCR(2+i);
        end
    end

    % plot the bars
    figure
    [h, order] = plotBarPlusDiff(plotdata, 'BaseValue', base, ...
        'colors', colors, 'newFigure', false);

    % add shaded regions
    num_groups = 2 + length(A.lambdaS); % inp, glp, dec times the number of lambdaS values
    len_group = length(A.lambdaC);
    for i = 1:num_groups
        for j = 1:length(A.input_SDRs)
            x1 = h(1 + (i-1)*len_group).XData(j) - h(1 + (i-1)*len_group).BarWidth/2;
            x2 = h(i*len_group).XData(j) + h(i*len_group).BarWidth/2;
            xregion(x1, x2,  'FaceColor', mapcolors(i, :), 'LineStyle', 'none', 'FaceAlpha', 0.2)
        end
    end
    
    % add legend etc.
    h = h(order(:));
    legend(h(1:length(legstr)), legstr, 'NumColumns', 1, 'Location', 'eastoutside')
    xlim([h(1).XData(1)-0.1, h(end).XData(end)+0.1])
    xticks(1:length(A.input_SDRs))
    xticklabels(A.input_SDRs)
    set(gca, 'ygrid', 'on')
    set(gca, 'yminorgrid', 'on')
    box('on')
    if metric == 1
        ylim([-1, max(plotdata(:))+1])
    else
        ylim([-4, 0])
    end
    xlabel('input SDR (dB)')
    
    % add title
    switch metric
        case 1
            title('$$\Delta$$SDR (dB), mean from all signals')
            ylabel('$$\Delta$$SDR (dB)')
        case 2
            title('PEMO-Q ODG, mean from all signals')
            ylabel('ODG')
        case 3
            title('PEAQ ODG, mean from all signals')
            ylabel('ODG')
    end

    % create a fake figure to use for matlab2tikz
    % ax = gca;
    % figure
    % b = bar(zeros(2, length(legstr)), 'grouped', 'FaceColor', 'flat');
    % for i = 1:length(b)
    %     b(i).CData = colors(i, :);
    % end
    % legend(legstr, 'location', 'eastoutside')
    % xlabel(ax.XLabel.String)
    % ylabel(ax.YLabel.String)
    % title(ax.Title.String)
    % set(gca, 'xlim', ax.XLim)
    % set(gca, 'ylim', ax.YLim)
    % set(gca, 'xtick', ax.XTick)
    % set(gca, 'xticklabel', ax.XTickLabel)
    % set(gca, 'ygrid', ax.YGrid)
    % set(gca, 'yminorgrid', ax.YMinorGrid)
    % set(gca, 'box', ax.Box)
end

%% make the figures docked
FigList = findobj(allchild(0), "flat", "Type", "figure");
FigList = flip(FigList);
for f = FigList
    set(f, 'WindowStyle', 'docked')
end

%% create color string for tikz
str = strings(size(colors, 1), 1);
for i = 1:size(colors, 1)
    str(i) = sprintf("\\definecolor{mycolor%d}{rgb}{%.5f,%.5f,%.5f}", ...
        i, colors(i, 1), colors(i, 2), colors(i, 3));
end