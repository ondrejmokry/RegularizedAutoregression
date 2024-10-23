% the script takes the EBU SQAM dataset (10 signals approximately 7 seconds
% long) and creates a shortened version for preliminary testing
%
% Date: 21/07/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

clear
clc
close all

%% load the EBU SQAM set
S = load('../survey toolbox/Sounds/Sound_database.mat');

% known sampling frequency
fs = 44100;

%% read the signal names
names = fields(S); % this returns a cell array of the names

%% crop the audio
% set the parameters
start = round(0.8*fs);
N = fs*0.6; % signal length
fi = round(N/20); % fade-in / fade-out length
cosinus = cos(linspace(0,pi/2,fi)').^2; % fade-in / fade-out shape

% crop and fade each signal
I = length(names); % number of signals
for i = 1:I
    data = S.(names{i});
    data = data(start+1:start+N);
    data = data/max(abs(data));
    data(1:fi) = data(1:fi).*(1-cosinus);
    data(end-fi+1:end) = data(end-fi+1:end).*cosinus;
    NS.(names{i}) = data;
end

%% plot the signals and spectrograms
figure
for i = 1:I
    % time domain
    subplot(5,4,2*i-1)
    plot(1000*(1:N)/fs,NS.(names{i}))
    title(names{i},'Interpreter','none')
    xlim(1000*[1 N]/fs)
    xlabel('Time (ms)')
    
    % spectrogram
    subplot(5,4,2*i)
    spectrogram(NS.(names{i}),[],[],[],fs,'yaxis')
    title(names{i},'Interpreter','none')
end

%% save
NS.fs = fs;
NS.names = names;
save('small_set.mat','-struct','NS')

%% sound
for i = 1:I
    fprintf('Playing sound: %s\n',names{i})
    pause(1)
    soundsc(NS.(names{i}),fs)
    pause(1)
end