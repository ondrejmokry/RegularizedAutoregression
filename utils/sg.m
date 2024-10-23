function sg(signal, fs, varargin)
% SG is a customized function for plotting the (power) spectrogram of the
% input signal, sampled at the input sampling frequency fs, with possible
% variable arguments:
%
% name        ('')     name to be displayed as a title of the figure
% solver      ('spt')  the toolbox used to compute the STFT, either
%                      Signal Processing Toolbox ('spt'), or LTFAT
%                      ('ltfat')
% reassigned  (true)   switch of the reassignment of the spectrogram
% w           (128)    window length in ms
% croplims    (true)   crops color axis
%
% Date: 03/10/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: xmokry12@vut.cz

% create the parser
pars = inputParser;

% add optional name-value pairs
addParameter(pars, 'name', '')
addParameter(pars, 'solver', 'spt')
addParameter(pars, 'reassigned', true)
addParameter(pars, 'w', 128)
addParameter(pars, 'croplims', true)

% parse
parse(pars, varargin{:})

% check signal size (should be a column vector)
if size(signal, 1) < size(signal, 2)
    signal = transpose(signal);
end

% compute and plot
w = floor(pars.Results.w*fs/1000);
M = w;
if strcmpi(pars.Results.solver, 'ltfat')
    % use LTFAT to compute the DGT
    if pars.Results.reassigned
        alpha = 0.75;
    else
        alpha = 0;
    end
    xres = 1600;
    yres = 1600;
    resgram(signal, fs, 'sharp', alpha, 'wlen', w, 'xres', xres, 'yres', yres)
else
    % use Signal Processing Toolbox to compute the DGT
    if pars.Results.reassigned
        [~, F, T, P] = spectrogram(signal, w, 3*w/4, M, fs, 'yaxis', 'reassigned');
    else
        [~, F, T, P] = spectrogram(signal, w, 3*w/4, M, fs, 'yaxis');
    end
    imagesc(T, F, 10*log10(P+eps))
    axis xy
end

% axis labels
xlabel('Time (s)')
ylabel('Frequency (Hz)')

% colormap
load('cmap', 'cmap')
colormap(cmap)

% colorbar
h = colorbar;
h.Label.String = 'Power/frequency (dB/Hz)';
h.Label.Interpreter = 'latex';

% plot title
if ~isempty(pars.Results.name)
    if contains(pars.Results.name, '_')
        title(pars.Results.name, 'Interpreter', 'none')
    else
        title(pars.Results.name, 'Interpreter', 'latex')
    end
end

% limits
ylim([0 min(fs/2, 12e3)])
if pars.Results.croplims
    lims = clim;
    clim(lims(1) + [0.6, 0.9]*(lims(2)-lims(1)))
end

end