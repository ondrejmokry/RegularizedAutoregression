function sg(signal, fs, options)
% SG is a customized function for plotting the (power) spectrogram of the
% input signal, sampled at the input sampling frequency fs, with possible
% variable arguments:

arguments
    signal
    fs
    options.name (1,:) char = ''            % name to be displayed as a title of the figure
    options.solver (1,:) char = 'spt'       % the toolbox used to compute the STFT, either Signal Processing Toolbox ('spt'), or LTFAT ('ltfat')
    options.reassigned (1,1) logical = true % switch of the reassignment of the spectrogram
    options.w (1,1) double = 128            % window length in ms
    options.croplims (1,1) logical = true   % crops color axis
end
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

% check signal size (should be a column vector)
if size(signal, 1) < size(signal, 2)
    signal = transpose(signal);
end

% compute and plot
w = floor(options.w*fs/1000);
M = w;
if strcmpi(options.solver, 'ltfat')
    % use LTFAT to compute the DGT
    if options.reassigned
        alpha = 0.75;
    else
        alpha = 0;
    end
    xres = 1600;
    yres = 1600;
    resgram(signal, fs, 'sharp', alpha, 'wlen', w, 'xres', xres, 'yres', yres)
else
    % use Signal Processing Toolbox to compute the DGT
    if options.reassigned
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
if ~isempty(options.name)
    if contains(options.name, '_')
        title(options.name, 'Interpreter', 'none')
    else
        title(options.name, 'Interpreter', 'latex')
    end
end

% limits
ylim([0 min(fs/2, 12e3)])
if options.croplims
    lims = clim;
    clim(lims(1) + [0.6, 0.9]*(lims(2)-lims(1)))
end

end