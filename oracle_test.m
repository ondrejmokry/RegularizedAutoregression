% the script tests the inpainting / declipping using Janssen algorithm or
% GLP and compares the progression of AR coefficients to the coefficients
% of the ground truth signal
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

clear
clc
close all
%#ok<*UNRCH>

addpath('utils')

%% settings
p        = 128;  % order of the AR process 
N        = 2048; % length of the signal
theta    = 0.2;  % threshold of the hard clipping
lambda   = [1, Inf]; % regularization parameters (AR coefficients, signal)
savedata = true;
simulate = true;

% construction of the matrix AA or XX in the subproblems
% ('toeplitz', 'xcorr', 'conv' and 'fft')
mat = 'toeplitz';

% used algos
methods = {'inpainting', 'inpainting (c)', 'janssen', 'janssen (c)', 'janssen (s)', 'janssen (c, s)', 'glp', 'glp (c)'};
fprintf('----------------------------\n')
fprintf('notation:                   \n')
fprintf(' (c) ... coef. acceleration \n')
fprintf(' (s) ... signal acceleration\n')
fprintf(' (c, s) ... both            \n')
fprintf('----------------------------\n')

%% load signal, normalize and crop it
% violin .... own sample
% piano ..... Beethoven Piano Sonata 21, 1st movement, bars 78-84
% quartet ... Beethoven Quartet Op. 18 No. 3, first movement, bars 156-162
% symphony .. Beethoven Symphony No. 9, finale opening bars
[signal, fs] = audioread('signals/violin.wav');

% convert to mono if necessary
signal = mean(signal, 2);

% shorten the signal
start  = round(0.95*fs);
signal = signal(start:start+N-1);

% in the simmulated case, we take the coefficients of the real signal and
% simulate an AR process based on these coefficients
if simulate
    c = lpc(signal, p);
    rng(0)
    noise = randn(length(signal)+p, 1);
    signal = filter(1, c, noise);
    signal = signal(p+1:end);
end

% normalize the signal
signal = signal/max(abs(signal));

%% create the degraded versions and the associated projection operators
masks.R  = abs(signal) < theta;
masks.U  = signal >= theta;
masks.L  = signal <= -theta;
degraded = max(-theta, min(signal, theta));

%% inpaint and declip the signal using the oracle coefficients
o_coef = getcoef(signal, p, lambda(1));
o_inpainted = inpaint(degraded, masks, o_coef);
o_declipped = declip(degraded, masks, o_coef);

%% initalize settings and results
settings = struct();
results = struct();

%% test everything
for variant = 1:5
    
    % if variant < 3
    %     DRmaxit = 'progressive';
    %     iterations = 10;
    % else
    %     DRmaxit = 1000;
    %     iterations = 5;
    % end

    DRmaxit = 1000;
    iterations = 100;

    switch variant
        case 1, coefextra = 0; sigextra = 0; linesearch = 0;
        case 2, coefextra = 0; sigextra = 0; linesearch = 1;
        case 3, coefextra = 1; sigextra = 0; linesearch = 0;
        case 4, coefextra = 0; sigextra = 1; linesearch = 0;
        case 5, coefextra = 1; sigextra = 1; linesearch = 0;
    end
    
    %% save the settings
    settings(variant).p          = p;
    settings(variant).lambda     = lambda;
    settings(variant).mat        = mat;
    settings(variant).DRmaxit    = DRmaxit;
    settings(variant).iterations = iterations;
    settings(variant).coefextra  = coefextra;
    settings(variant).sigextra   = sigextra;
    settings(variant).linesearch = linesearch;
       
    fprintf('============================\n')
    fprintf('coefficient extrapolation: %d\n', coefextra)
    fprintf('signal extrapolation:      %d\n', sigextra)
    fprintf('linesearch:                %d\n', linesearch)
    fprintf('============================\n')
    
    %% initialization
    signals = NaN(N, iterations, length(methods));
    obj = NaN(iterations, length(methods));
    time = NaN(iterations, length(methods));

    %% set up the profiler
    % profile clear
    % profile -memory on

    %% inpaint
    fprintf('Inpainting: ')
    tic
    [signals(:, :, 1), obj(:, 1), time(:, 1)] = janssen('inpainting', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit);
    toc

    fprintf('Inpainting (c): ')
    tic
    [signals(:, :, 2), obj(:, 2), time(:, 2)] = janssen('inpainting', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit, ...
        'coefaccel', true);
    toc

    %% declip using Janssen
    fprintf('Janssen: ')
    tic
    [signals(:, :, 3), obj(:, 3), time(:, 3)] = janssen('declipping', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit);
    toc

    fprintf('Janssen (c): ')
    tic
    [signals(:, :, 4), obj(:, 4), time(:, 4)] = janssen('declipping', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit, ...
        'coefaccel', true);
    toc
    
    fprintf('Janssen (s): ')
    tic
    [signals(:, :, 5), obj(:, 5), time(:, 5)] = janssen('declipping', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit, ...
        'sigaccel', true);
    toc

    fprintf('Janssen (c, s): ')
    tic
    [signals(:, :, 6), obj(:, 6), time(:, 6)] = janssen('declipping', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit, ...
        'sigaccel', true, ...
        'coefaccel', true);
    toc

    %% declip using GLP
    fprintf('GLP: ')
    tic
    [signals(:, :, 7), obj(:, 7), time(:, 7)] = janssen('glp', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit);
    toc

    fprintf('GLP (c): ')
    tic
    [signals(:, :, 8), obj(:, 8), time(:, 8)] = janssen('glp', degraded, masks, lambda, p, iterations, ...
        'coefextra', coefextra, ...
        'sigextra', sigextra, ...
        'linesearch', linesearch, ...
        'saveall', true, ...
        'DRmaxit', DRmaxit, ...
        'coefaccel', true);
    toc

    %% view the profiler
    % profile viewer
    
    %% compute coefficients during iterations (yes, I know it could be output of the algorithms...)
    coefs     = NaN(p+1, iterations, length(methods));
    relatives = NaN(iterations-1, length(methods));
    distances = NaN(iterations, length(methods));
    SDRs      = NaN(iterations, length(methods));
    for method = 1:length(methods)
        for i = 1:iterations
            % the following condition ensures that the iteration was really
            % computed during the frame procedure (i.e., it was not stopped by
            % the PSD check)
            if ~isnan(time(i, method))
                coefs(:, i, method) = getcoef(signals(:, i, method), p, lambda(1));
                if i > 1
                    relatives(i-1, method) = norm(coefs(:, i, method)-coefs(:, i-1, method))/norm(coefs(:, i, method));
                end
                distances(i, method) = norm(coefs(:, i, method)-o_coef);
                SDRs(i, method) = snr(signal(~masks.R), signal(~masks.R)-signals(~masks.R, i, method));
            end
        end
    end
    
    %% save the solution
    results(variant).coefs     = coefs;
    results(variant).distances = distances;
    results(variant).obj       = obj;
    results(variant).relatives = relatives;
    results(variant).SDRs      = SDRs;
    results(variant).signals   = signals;
    results(variant).time      = time;
    
end

%% save the data
% save the current time
d = datetime('now');

% save using automatically defined name
if savedata
    save(fname, ...
        'degraded', ...
        'fs', ...
        'masks', ...
        'methods', ...
        'N', ...
        'o_coef', ...
        'o_inpainted', ...
        'o_declipped', ...  
        'results', ...
        'settings', ...
        'signal', ...
        'theta')
end

%% plot
for variant = 1:length(settings)
    
    % solutions
    figure('visible', 'off')
    C = colororder;
    colororder([C; 0.6*C(1, [2 3 1])])
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

    % title
    sgtitle(sprintf('coefficient extrapolation: %d\nsignal extrapolation: %d\nlinesearch: %d', ...
        settings(variant).coefextra, ...
        settings(variant).sigextra, ...
        settings(variant).linesearch))
    
    % make the figure full-screen
    set(gcf, 'WindowState', 'maximized', 'visible', 'on')
    
end

%% functions
function [coef, objective] = getcoef(signal, p, lambda)

    X  = toeplitz([signal; zeros(p, 1)], [signal(1), zeros(1, p)]);
    XX = X'*X;

    % Douglas-Rachford
    DR.lambda = 1;
    DR.gamma  = 10;
    DR.y0     = lpc(signal, p)';
    DR.maxit  = 1000;
    DR.tol    = -Inf;
    DR.dim    = p+1;

    % the function g
    DR.g      = @(a) lambda*norm(a, 1);
    soft      = @(a, t) sign(a).*max(abs(a)-t, 0);
    DR.prox_g = @(a, t) [1; soft(a(2:end), lambda*t)]; 

    % the function f
    DR.f      = @(a) 0.5*norm(X*a)^2;
    invmat    = inv(eye(p+1) + DR.gamma*XX);
    DR.prox_f = @(a, t) invmat*a; %#ok<MINV>

    % solution
    if nargout > 1
        [coef, objective] = DouglasRachford(DR, DR);
    else
        coef = DouglasRachford(DR, DR);
    end
    
end

function inpainted = inpaint(degraded, masks, coef)

    % prepare the dimensions
    N         = length(degraded);
    p         = length(coef)-1;
    
    % perform inpainting step as in Audio Inpainting Toolbox
    inpainted = degraded;
    indmiss   = find(~masks.R);
    indobs    = find(masks.R);
    IAA       = abs(repmat(indmiss, 1, N)-repmat(1:N, length(indmiss), 1));
    IAA1      = IAA <= p;
    AA        = zeros(size(IAA));
    b         = coef'*hankel(coef', [coef(end), zeros(1, p)]);
    AA(IAA1)  = b(IAA(IAA1)+1);
    R         = chol(AA(:, indmiss));
    inpainted(~masks.R) = -R\(R'\(AA(:, indobs)*degraded(indobs)));
    
end

function declipped = declip(clipped, masks, coef)

    N = length(clipped);
    proj = @(x) clipped .* masks.R + ...
            max(x, clipped) .* masks.U + ...
            min(x, clipped) .* masks.L;

    % prepare the matrix A
    A  = toeplitz([coef; zeros(N-1, 1)], [coef(1), zeros(1, N-1)]);
    AA = A'*A;
    
    % set the parameters of the algorithm
    DR.lambda = 1;
    DR.gamma  = 10;
    DR.y0     = clipped;
    DR.maxit  = 2000;
    DR.tol    = -Inf;

    % precompute the inversion
    invmat = inv(eye(N) + DR.gamma*AA);

    % set the parameters of the model
    DR.prox_f = @(x, t) invmat*x; %#ok<MINV>
    DR.prox_g = @(x, t) proj(x); 
    DR.dim    = N;

    % solve the task
    declipped = DouglasRachford(DR, DR);

end