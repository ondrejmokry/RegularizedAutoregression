clear
clc
close all

load("demo_times.mat")

methods = fieldnames(methodsTable);

%% compute time per iteration
for m = 1:length(methods)
    timeperit = zeros(height(resultsTable), 1);
    for i = 1:height(resultsTable)
        numit = sum(~isnan(resultsTable.time{i}(:, m)));
        timeperit(i) = resultsTable.time{i}(numit, m)/numit; % number of iterations
    end
    methodsTable.(methods{m}) = addvars(methodsTable.(methods{m}), timeperit, 'NewVariableNames', 'timeperit');
end

%% convert data types
types = [ ...
    "double", "double", "double", "double", "double", "double", ... % first set of parameters
    "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", ... % flags
    "double", "double", "double", ... % second set of parameters
    "cell", "cell", "cell", "cell" ... % results
    ];
for i = 1:18
    resultsTable.(i) = cell2mat(resultsTable{:, i});
end
resultsTable.Properties.VariableTypes = types;
for m = 1:length(methods)
    for i = 1:22
        methodsTable.(methods{m}).(i) = cell2mat(methodsTable.(methods{m}){:, i});
    end
    methodsTable.(methods{m}).Properties.VariableTypes = ...
        [types(1:18), "double", "double", "double", "double", "double"];
end

%% plot settings
methodnames = strrep(methods, '_', ' ');
suffixes = {'', ' sparse'};
colors = lines(4);
styles = {'-'; '--'};
widths = [0.5, 1.0];

figure
tiledlayout(1, 3)

%% dependence on model order (p) for fixed N = 4096
fixN = 4096;
fixtheta = 0.2;
rows = [ ...
    (resultsTable.N == fixN) .* (resultsTable.theta == fixtheta) .* (resultsTable.lambdaC == 0), ...
    (resultsTable.N == fixN) .* (resultsTable.theta == fixtheta) .* (resultsTable.lambdaC > 0) ...
    ];
rows = logical(rows);

nexttile
for m = 1:length(methods)
    for i = 1:2
        x = methodsTable.(methods{m}).p(rows(:, i));
        y = methodsTable.(methods{m}).timeperit(rows(:, i));
        plot(x, y, ...
            Color=colors(m, :), ...
            LineStyle=styles{i}, ...
            LineWidth=widths(i), ...
            DisplayName=[methodnames{m}, suffixes{i}])
        hold on
    end
end
grid on
grid minor
xlabel("model order p")
ylabel("elapsed time per outer iteration (s)")
title("dependence on model order (p)", "for fixed N = 4096")

%% dependence on segment length (N) for fixed p = 256
fixp = 256;
fixtheta = 0.2;
rows = [...
    (resultsTable.p == fixp) .* (resultsTable.theta == fixtheta) .* (resultsTable.lambdaC == 0), ...
    (resultsTable.p == fixp) .* (resultsTable.theta == fixtheta) .* (resultsTable.lambdaC > 0), ...
    ];
rows = logical(rows);

nexttile
for m = 1:length(methods)
    for i = 1:2
        x = methodsTable.(methods{m}).N(rows(:, i));
        y = methodsTable.(methods{m}).timeperit(rows(:, i));
        plot(x, y, ...
            Color=colors(m, :), ...
            LineStyle=styles{i}, ...
            LineWidth=widths(i), ...
            DisplayName=[methodnames{m}, suffixes{i}])
        hold on
    end
end
grid on
grid minor
xlabel("segment length N")
ylabel("elapsed time per outer iteration (s)")
title("dependence on segment length (N)", "for fixed p = 256")

%% dependence on clipping threshold (theta) for fixed N = 4096 and p = 256
fixN = 4096;
fixp = 256;
rows = [ ...
    (resultsTable.N == fixN) .* (resultsTable.p == fixp) .* (resultsTable.lambdaC == 0), ...
    (resultsTable.N == fixN) .* (resultsTable.p == fixp) .* (resultsTable.lambdaC > 0), ...
    ];
rows = logical(rows);

nexttile
for m = 1:length(methods)
    for i = 1:2
        x = methodsTable.(methods{m}).theta(rows(:, i));
        y = methodsTable.(methods{m}).timeperit(rows(:, i));
        plot(x, y, ...
            Color=colors(m, :), ...
            LineStyle=styles{i}, ...
            LineWidth=widths(i), ...
            DisplayName=[methodnames{m}, suffixes{i}])
        hold on
    end
end
lgd = legend;
lgd.Layout.Tile = "east";
grid on
grid minor
xlabel("clipping threshold theta")
ylabel("elapsed time per outer iteration (s)")
title("dependence on clipping threshold (theta)", "for fixed N = 4096 and p = 256")
