% the script plots the data generated by oracle_test.m

clear
clc
close all
%#ok<*UNRCH>
addpath('../utils')

%% load the file
load('oracle_test_01.mat')

%% setup
% methods:
% 1: inpainting
% 2: inpainting (c)
% 3: janssen
% 4: janssen (c)
% 5: janssen (s)
% 6: janssen (c, s)
% 7: glp
% 8: glp (c)
plotmethods = [1, 3, 7];

% possibility to crop last sample from GLP metrics because no rectification
% is done in the last iteration and this might affect the results a lot
cropglp = true;

%% choose data based on setup
% only some methods
for variant = 1:length(settings)
    results(variant).coefs = results(variant).coefs(:, :, plotmethods); %#ok<*SAGROW>
    results(variant).distances = results(variant).distances(:, plotmethods);
    results(variant).obj = results(variant).obj(:, plotmethods);
    results(variant).SDRs = results(variant).SDRs(:, plotmethods);
    results(variant).relatives = results(variant).relatives(:, plotmethods);
    results(variant).signals = results(variant).signals(:, :, plotmethods);
    results(variant).time = results(variant).time(:, plotmethods);
end
methods = methods(plotmethods);

% crop GLP data
if cropglp
    for m = 1:length(methods)
        if contains(methods{m}, 'glp')
            for variant = 1:length(settings)
                results(variant).distances(end, m) = NaN;
                results(variant).obj(end, m) = NaN;
                results(variant).relatives(end, m) = NaN;
                results(variant).SDRs(end, m) = NaN;
                results(variant).time(end, m) = NaN;
            end
        end
    end
end

%% plot
for variant = 1:length(settings)
    
    % solutions
    figure('visible', 'off')
    C = colororder;
    colororder(repmat(C(1:length(methods), :), 2, 1))
    subplot(2, 4, 1:4)
    h1 = plot(squeeze(results(variant).signals(:, end, :)));
    hold on
    if variant > 1
        plot(squeeze(results(1).signals(:, end, :)), ':')
    end
    h2 = plot(signal, 'color', 0.75*[1, 1, 1]);
    h3 = plot(o_inpainted, 'color', 0.50*[1, 1, 1]);
    h4 = plot(o_declipped, 'color', 0.25*[1, 1, 1]);
    legend([h1; h2; h3; h4], [methods, {'original', 'inpainted from original coefs', 'declipped from original coefs'}])
    xlim([1 N])

    % objectives
    subplot(2, 4, 5)
    h5 = loglog(results(variant).time, results(variant).obj, '-x');
    if variant > 1
        hold on
        loglog(results(1).time, results(1).obj, ':')
    end
    xlabel('time (s)')
    ylabel('objective')
    legend(h5, methods)
    grid on
    grid minor
    ylim('padded')

    % SDRs
    subplot(2, 4, 6)
    h6 = loglog(results(variant).time, results(variant).SDRs, '-x');
    if variant > 1
        hold on
        loglog(results(1).time, results(1).SDRs, ':')
    end
    xlabel('time (s)')
    ylabel('SDR (dB)')
    legend(h6, methods, 'location', 'southeast')
    grid on
    grid minor
    ylim('padded')

    % coef distances
    subplot(2, 4, 7)
    h7 = loglog(results(variant).time, results(variant).distances, '-x');
    if variant > 1
        hold on
        loglog(results(1).time, results(1).distances, ':')
    end
    xlabel('time (s)')
    ylabel('norm of coef difference')
    legend(h7, methods)
    grid on
    grid minor
    ylim('padded')

    % plot coef relative distances
    subplot(2, 4, 8)
    h8 = loglog(results(variant).time(2:end, :), results(variant).relatives, '-x');
    if variant > 1
        hold on
        loglog(results(1).time(2:end, :), results(1).relatives, ':')
    end
    xlabel('time (s)')
    ylabel('relative norms of coefficients during iterations')
    legend(h8, methods)
    grid on
    grid minor
    ylim('padded')

    % title
    sgtitle(sprintf('coefficient extrapolation: %d\nsignal extrapolation: %d\nlinesearch: %d', ...
        settings(variant).coefextra, ...
        settings(variant).sigextra, ...
        settings(variant).linesearch))
    
    % make the figure full-screen
    set(gcf, 'WindowState', 'maximized', 'visible', 'on')
    
end