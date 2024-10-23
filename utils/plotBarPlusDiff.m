function [h, order] = plotBarPlusDiff(stackData, varargin)
% Plot a set of stacked bars and group them according to labels provided.
% Based on the function plotBarStackGroups from
% https://www.mathworks.com/matlabcentral/fileexchange/32884-plot-groups-of-stacked-bars
%
% This version is meant to display a comparison between two sets of grouped
% data.
%
% Params: 
%    stackData is a 3D array organized such that
%       size(stackData, 1) = num of groups
%       size(stackData, 2) = num of stacks per group
%       size(stackData, 3) = 2 = num of elements per stack
%    groupLabels, stackLabels and elementLabels are cells of strings
%       defining the corresponding labels for axis label or legend
%    newFigure defines whether the plot will be in a new figure
%
% Copyright 2011 Evan Bollig (bollig@scs.fsu.edu)
% Edited 2021 Ondrej Mokry (ondrej.mokry@mensa.cz)

% read the size of the input
numGroups   = size(stackData, 1);
numStacks   = size(stackData, 2);
numElements = size(stackData, 3);

% create the parser
pars = inputParser;

% add optional name-value pairs
addParameter(pars, 'groupLabels', [])
addParameter(pars, 'stackLabels', [])
addParameter(pars, 'elementLabels', [])
addParameter(pars, 'newFigure', true)
addParameter(pars, 'legend', true)
addParameter(pars, 'colors', 'default')
addParameter(pars, 'baseValue', 0)

% parse
parse(pars, varargin{:})
groupLabels = pars.Results.groupLabels;
stackLabels = pars.Results.stackLabels;
if isempty(stackLabels)
    stackLabels = cell(numStacks, 1);
end
elementLabels = pars.Results.elementLabels;
if isempty(elementLabels)
    elementLabels = cell(numElements, 1);
end

% decide on the colors
if isfloat(pars.Results.colors)
    colors = pars.Results.colors;
else
    colors = parula(numStacks);
end

% set up the bins
groupBins = 1:numGroups;
maxGroupWidth = 0.75;

if pars.Results.newFigure
    figure
end
hold on
permsFirstGroup = repmat(1:numElements, numStacks, 1);
for i = 1:numStacks
    Y = squeeze(stackData(:, i, :));
    perms = repmat(1:numElements, numGroups, 1);
    
    % compute the difference and swap the data if necessary
    for j = 1:numGroups
        if Y(j, 2) >= 0
            if Y(j, 1) >= 0
                [ Y(j, :), perms(j, :) ] = sort(Y(j, :));
                % plot the smaller value Y(j, 1) and then the difference
                Y(j, 2) = Y(j, 2) - Y(j, 1);
            end
        else
            if Y(j, 1) < 0
                if pars.Results.baseValue >= 0
                    [ Y(j, :), perms(j, :) ] = sort(Y(j, :), 'descend');
                    % plot the bigger value Y(j, 1) and then the difference
                    Y(j, 2) = Y(j, 2) - Y(j, 1);
                else
                    [ Y(j, :), perms(j, :) ] = sort(Y(j, :), 'descend');
                    % plot the smaller value Y(j, 1) and then the
                    % difference, but switch the colors
                    Y(j, 2) = Y(j, 2) - Y(j, 1);
                    perms(j, :) = flip(perms(j, :));
                end
            end
        end
        
        % in other cases, one value is >= 0 and the second one is < 0, so
        % both are plotted in their original order without computing the
        % difference
    end
    
    % save the permutations of the first group in order to correctly
    % assign the legend
    permsFirstGroup(i, :) = perms(1, :);
    
    % center the bars
    internalPosCount = i-((numStacks+1)/2);
    
    % offset the group draw positions:
    groupDrawPos = (internalPosCount)*(6/5*maxGroupWidth/numStacks) + groupBins;
    
    h(i, :) = bar(Y, 'stacked', 'BaseValue', pars.Results.baseValue); %#ok<AGROW>
    set(h(i, :), 'BarWidth', maxGroupWidth/numStacks);
    set(h(i, :), 'XData', groupDrawPos);
    for j = 1:numElements
        h(i, j).FaceColor = 'flat';
    end
    for j = 1:numGroups
        for k = 1:numElements
            mixture = (numElements - k + 1)/numElements;
            mixture = mixture^2;
            h(i, perms(j, k)).CData(j, :) = ...
                colors(i, :) * mixture + [1 1 1] * (1-mixture);
        end
    end
end
hold off
set(gca, 'XTickMode', 'manual')
set(gca, 'XTick', 1:numGroups)
if ~isempty(groupLabels)
    set(gca, 'XTickLabelMode', 'manual')
    set(gca, 'XTickLabel', groupLabels)
end

% define legend
legcell = cell(numStacks, numElements);
for i = 1:numStacks
    for j = 1:numElements
        legcell{i, j} = [char(stackLabels{i}), char(elementLabels{j})];
    end
end

% reorder the legend according to the first group plotted
order = [(1:numStacks)', (numStacks+1:2*numStacks)'];
for i = 1:numStacks
    order(i, :) = order(i, permsFirstGroup(i, :));
end

% add the legend to the plot
if pars.Results.legend
    legend(h(order(:)), legcell, 'NumColumns', numElements)
end

end