% the script plots the data generated by the script 'iteration_tradeoff.m'
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

clear
clc
close all

addpath('../utils')

%% load the data
load('iteration_tradeoff_01.mat')
chosenlamb = [1, 4];

%% prepare the grid
[X, Y] = meshgrid(1:maxit, DRmaxit);

%% set parameters and initialize
I = length(chosenlamb);
levels = 40;
[A, B] = sbplts(I);
f1 = figure('visible', 'off');
f2 = figure('visible', 'off');
f3 = figure('visible', 'off');
f4 = figure('visible', 'off');

% plotfn = @(x, y, z) contourf(x, y, z, levels, 'EdgeAlpha', 0);
plotfn = @(x, y, z) contour(x, y, z, levels);

%% plot
for i = 1:I
    % SDR obtained via declipFrameJanssen
    set(0, 'CurrentFigure', f1)
    subplot(A, B, i)
    plotfn(X, Y, squeeze(SDRone(chosenlamb(i), :, :))')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('outer iterations')
    ylabel('DR iterations')
    title(['lambda = ', num2str(lambda(chosenlamb(i)))])
    colorbar
    
    % objective obtained via declipFrameJanssen
    set(0, 'CurrentFigure', f2)
    subplot(A, B, i)
    plotfn(X, Y, log10(squeeze(OBJone(chosenlamb(i), :, :)))')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('outer iterations')
    ylabel('DR iterations')
    title(['lambda = ', num2str(lambda(chosenlamb(i)))])
    colorbar
    
    % SDR obtained via declipFrameGLP
    set(0, 'CurrentFigure', f3)
    subplot(A, B, i)
    plotfn(X, Y, squeeze(SDRtwo(chosenlamb(i), :, :))')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('outer iterations')
    ylabel('DR iterations')
    title(['lambda = ', num2str(lambda(chosenlamb(i)))])
    colorbar
    
    % objective obtained via declipFrameGLP
    set(0, 'CurrentFigure', f4)
    subplot(A, B, i)
    plotfn(X, Y, log10(squeeze(OBJtwo(chosenlamb(i), :, :)))')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('outer iterations')
    ylabel('DR iterations')
    title(['lambda = ', num2str(lambda(chosenlamb(i)))])
    colorbar    
end

%% add figure titles
sgtitle(f1, {'Janssen, SDR', ...
    sprintf('gammaC = %.1f, gammaS = %.1f', gammaC, gammaS)})
sgtitle(f2, {'Janssen, log10 of objective', ...
    sprintf('gammaC = %.1f, gammaS = %.1f', gammaC, gammaS)})
sgtitle(f3, {'GLP, SDR', sprintf('gammaC = %.1f', gammaC)})
sgtitle(f4, {'GLP, log10 of objective', sprintf('gammaC = %.1f', gammaC)})

%% make figures maximized and visible
for f = [f1, f2, f3, f4]
    % set(f, 'WindowState', 'maximized', 'visible', 'on')
    set(f, 'visible', 'on')
end