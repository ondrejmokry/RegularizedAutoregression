clear
clc
close all
%#ok<*UNRCH>

addpath("utils")

savetimes = false;

%% settings
p            = 256;   % order of the AR process 
N            = 4096;  % length of the signal
theta        = 0.8;   % threshold of the hard clipping
iterations   = 10;    % Janssen iterations
lambdaC      = 0.1;   % parameter of AR coefficient regularization
lambdaS      = 10;    % parameter of incosistent signal regularization
simulate     = true;  % use simulated signal?
plottime     = false; % use time as horizontal axis in plots
segmentation = false; % process the signal segment-wise

% some more technical settings
coefextra  = false;
sigextra   = false;
coefaccel  = true;
sigaccel   = true;
linesearch = false;
plotLS     = false;
gammaC     = 0.1;
gammaS     = 10;
DRmaxit    = 1000;

% construction of the matrix AA or XX in the subproblems
% ("toeplitz", "xcorr", "conv" and "fft")
mat = "toeplitz";

% used algos
methods = ["inpainting", "GLP", "inconsistent declipping", "consistent declipping"];

%% load signal, normalize and crop it
% violin .... own sample
% piano ..... Beethoven Piano Sonata 21, 1st movement, bars 78-84
% quartet ... Beethoven Quartet Op. 18 No. 3, first movement, bars 156-162
% symphony .. Beethoven Symphony No. 9, finale opening bars
[signal, fs] = audioread("signals/violin.wav");

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

%% create the degraded version
masks.R  = abs(signal) < theta;
masks.U  = signal >= theta;
masks.L  = signal <= -theta;
degraded = max(-theta, min(signal, theta));
    
%% initialization
signals = NaN(N, iterations, length(methods));
obj = NaN(iterations, length(methods));
time = NaN(iterations, length(methods));

%% processing
% inpainting
fprintf("Inpainting: ")
[signals(:, :, 1), O, T] = janssen("inpainting", ...
    degraded, masks, lambdaC, p, iterations, ...
    "segmentation", segmentation, ...
    "DRmaxit", DRmaxit, "decompose", true, "mat", mat, ...
    "gammaC", gammaC, "gammaS", gammaS, ...
    "coefaccel", coefaccel, "sigaccel", sigaccel, ...
    "coefextra", coefextra, "sigextra", sigextra, ...
    "linesearch", linesearch, "plotLS", plotLS, "saveall", true);
obj(:, 1) = sum(O, 2);
time(:, 1) = sum(T, 2);
fprintf("elapsed time %.1f seconds\n", time(end, 1))

% declipping using GLP
fprintf("GLP: ")
[signals(:, :, 2), O, T] = janssen("glp", ...
    degraded, masks, lambdaC, p, iterations, ...
    "segmentation", segmentation, ...
    "DRmaxit", DRmaxit, "decompose", true, "mat", mat, ...
    "gammaC", gammaC, "gammaS", gammaS, ...
    "coefaccel", coefaccel, "sigaccel", sigaccel, ...
    "coefextra", coefextra, "sigextra", sigextra, ...
    "linesearch", linesearch, "plotLS", plotLS, ...
    "saveall", true, "verbose", false);
obj(:, 2) = sum(O, 2);
time(:, 2) = sum(T, 2);
fprintf("elapsed time %.1f seconds\n", time(end, 2))

% declipping using consistent Janssen
fprintf("Inconsistent Janssen: ")
[signals(:, :, 3), O, T] = janssen("declipping", ...
    degraded, masks, [lambdaC, lambdaS], p, iterations, ...
    "segmentation", segmentation, ...
    "DRmaxit", DRmaxit, "decompose", true, "mat", mat, ...
    "gammaC", gammaC, "gammaS", gammaS, ...
    "coefaccel", coefaccel, "sigaccel", sigaccel, ...
    "coefextra", coefextra, "sigextra", sigextra, ...
    "linesearch", linesearch, "plotLS", plotLS, ...
    "saveall", true, "verbose", false);
obj(:, 3) = sum(O, 2);
time(:, 3) = sum(T, 2);
fprintf("elapsed time %.1f seconds\n", time(end, 3))

% declipping using consistent Janssen
fprintf("Consistent Janssen: ")
[signals(:, :, 4), O, T] = janssen("declipping", ...
    degraded, masks, [lambdaC, Inf], p, iterations, ...
    "segmentation", segmentation, ...
    "DRmaxit", DRmaxit, "decompose", true, "mat", mat, ...
    "gammaC", gammaC, "gammaS", gammaS, ...
    "coefaccel", coefaccel, "sigaccel", sigaccel, ...
    "coefextra", coefextra, "sigextra", sigextra, ...
    "linesearch", linesearch, "plotLS", plotLS, ...
    "saveall", true, "verbose", false);
obj(:, 4) = sum(O, 2);
time(:, 4) = sum(T, 2);
fprintf("elapsed time %.1f seconds\n", time(end, 4))

%% compute metrics
relatives = NaN(iterations-1, length(methods));
SDRs      = NaN(iterations, length(methods));
for method = 1:length(methods)
    for i = 1:iterations
        % the following condition ensures that the iteration was really
        % computed during the frame procedure (i.e., it was not stopped by
        % the PSD check)
        if ~isnan(time(i, method))
            if i > 1
                relatives(i-1, method) = norm(signals(:, i, method)-signals(:, i-1, method))/norm(signals(:, i, method));
            end
            SDRs(i, method) = snr(signal(~masks.R), signal(~masks.R)-signals(~masks.R, i, method));
        end
    end
end

%% save elapsed times
resultsNames = {'p', 'N', 'theta', 'iterations', 'lambdaC', 'lambdaS', ...
                 'simulate', 'plottime', 'segmentation', 'coefextra', ...
                 'sigextra', 'coefaccel', 'sigaccel', 'linesearch', ...
                 'plotLS', 'gammaC', 'gammaS', 'DRmaxit', ...
                 'time', 'obj', 'SDRs', 'relatives'};
resultsValues = {p, N, theta, iterations, lambdaC, lambdaS, ...
                  simulate, plottime, segmentation, coefextra, ...
                  sigextra, coefaccel, sigaccel, linesearch, ...
                  plotLS, gammaC, gammaS, DRmaxit, ...
                  time, obj, SDRs, relatives};
resultsTable = array2table(resultsValues, 'VariableNames', resultsNames);
methodsTable = struct;
for method = 1:length(methods)
    modname = strrep(methods(method), ' ', '_');
    methodValues = {p, N, theta, iterations, lambdaC, lambdaS, ...
                  simulate, plottime, segmentation, coefextra, ...
                  sigextra, coefaccel, sigaccel, linesearch, ...
                  plotLS, gammaC, gammaS, DRmaxit, ...
                  time(end, method), obj(end, method), SDRs(end, method), relatives(end, method)};
    methodsTable.(modname) = array2table(methodValues, 'VariableNames', resultsNames);
end
if isfile("results/demo_times.mat")
    D = load("results/demo_times.mat");
    D.resultsTable = sortrows([D.resultsTable; resultsTable], 1:18);
    for method = 1:length(methods)
        modname = strrep(methods(method), ' ', '_');
        D.methodsTable.(modname) = sortrows([D.methodsTable.(modname); methodsTable.(modname)], 1:18);
    end
    if savetimes
        save("results/demo_times.mat", "-struct", "D")
    end
elseif savetiemes
    save("results/demo_times.mat", "resultsTable", "methodsTable")
end


%% plot
figure
tiledlayout(2, 3)

% ignore some last values for GLP due to the rectification
obj(end, 2) = NaN;
relatives(end, 2) = NaN;

% deside on the horizontal axis
if plottime
    xaxis = time;
    xstr = "time (s)";
else
    xaxis = repmat((1:iterations)', [1, length(methods)]);
    xstr = "iteration";
end

% solutions
nexttile([1, 3])
h1 = plot(squeeze(signals(:, end, :)));
hold on
h2 = plot(degraded, "color", 0.75*[1, 1, 1]);
h3 = plot(signal, "color", 0.25*[1, 1, 1]);
legend([h3; h2; h1], ["original", "degraded", methods])
xlim([1 N])

% objectives
nexttile
plot(xaxis, obj, "-x")
xlabel(xstr)
ylabel("objective")
legend(methods)
grid on
grid minor
set(gca, "YScale", "log")
if plottime
    set(gca, "XScale", "log")
end

% SDRs
nexttile
plot(xaxis, SDRs, "-x")
xlabel(xstr)
ylabel("SDR (dB)")
legend(methods, "location", "southeast")
grid on
grid minor
if plottime
    set(gca, "XScale", "log")
end

% relative changes of solution
nexttile
plot(xaxis(2:end, :), relatives, "-x")
xlabel(xstr)
ylabel("relative solution change")
legend(methods)
grid on
grid minor
set(gca, "YScale", "log")
if plottime
    set(gca, "XScale", "log")
end
