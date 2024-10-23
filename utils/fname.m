function name = fname()
% fname sets the filename for saving the results based on the name of the
% script fname is called in
%
% it returns string formed from the script name and a numeric identifier
% such that no results are overwritten
%
% name also includes the folder path 'results/' as a prefix
%
% Date: 08/07/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% get the function call stack
stack = dbstack;

% read the second one (the first one is fname) and put a prefix
name = ['results/', stack(2).name, '_01.mat'];

% iterate to determine the numeric identifier
i = 1;
while isfile(['./', name])
    i = i + 1;
    name = [name(1:end-6), num2str(i,'%02d'), '.mat'];
end

