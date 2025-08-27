function [restored, objective, times] = janssen(method, signal, masks, lambda, p, maxit, options)
% janssen is a function that wraps four closely related reconstruction
% algorithms for audio declipping/inpainting:
% (1) inpainting using the Janssen algorithm [1] with sparse regularization
%     of the AR coefficients
% (2) declipping using the Janssen algorithm [1] with sparse regularization
%     of the AR coefficients and consistency constraints on the signal
% (3) declipping using the generalized linear prediction [2] with sparse
%     regularization of the AR coefficients
% (4) dequantization using the Janssen algorithm [1] with sparse 
%     regularization of the AR coefficients and consistency constraints on 
%     the signal
%
% implementation of the signal estimation step in (1) and (3) is taken from
% the Audio Inpainting Toolbox as described in [3]
%
% the subproblems of estimating the AR coefficients in (1)–(4) and
% the signal in (2) and (4) are solved using the Douglas–Rachford
% algorithm [4], which is possibly accelerated as described in [5]
%
% [1] A. Janssen, R. Veldhuis and L. Vries, "Adaptive interpolation of
%     discrete-time signals that can be modeled as autoregressive
%     processes," in IEEE Transactions on Acoustics, Speech, and Signal
%     Processing, vol. 34, no. 2, pp. 317-330, 1986, doi:
%     10.1109/TASSP.1986.1164824.
% [2] L. Atlas and C. P. Clark, "Clipped-waveform repair in acoustic
%     signals using generalized linear prediction," US patent US8126578B2, 
%     2007.
% [3] A. Adler, V. Emiya, M. G. Jafari, M. Elad, R. Gribonval and M. D.
%     Plumbley, "Audio Inpainting, " in IEEE Transactions on Audio, Speech, 
%     and Language Processing, vol. 20, no. 3, pp. 922-932, 2012, doi:
%     10.1109/TASL.2011.2168211.
% [4] P. Combettes and J.-C. Pesquet, "Proximal Splitting Methods in Signal
%     Processing," in Fixed-Point Algorithms for Inverse Problems in
%     Science and Engineering, vol. 49, pp. 185-212, 2009, doi:
%     10.1007/978-1-4419-9569-8_10.
% [5] I. Bayram, "Proximal Mappings Involving Almost Structured Matrices,"
%     in IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2264-2268, 
%     2015, doi: 10.1109/LSP.2015.2476381.
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%   objective     values of the objective function during iterations
%   times         cumulative computation time during iterations
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

arguments
    % Switch between the algorithms (1)–(3):
    % (1) "inpainting"
    % (2) "declipping"
    % (3) "glp"
    % (4) "dequantization"
    method (1,1) string {mustBeMember(method, ["inpainting", "declipping", "glp", "dequantization"])}
    
    % The input (degraded) signal
    signal (:,1) double
    
    % Masks:
    % masks.R ... mask of the reliable samples
    % masks.U ... mask of the samples clipped on the upper clipping level (needed only for (2), (3))
    % masks.L ... mask of the samples clipped on the lower clipping level (needed only for (2), (3))
    masks struct
    
    % Regularization parameters:
    % lambda(1) ... for the AR model estimation, if lambda(1) = 0, no regularization is used
    % lambda(2) ... for the signal estimation, if lambda(2) = Inf or length(lambda) = 1,
    %               the indicator function is used instead of the distance function
    lambda (1,:) double
    
    % Order of the AR model
    p (1,1) double {mustBePositive}
    
    % Number of iterations of the whole Janssen algorithm
    maxit (1,1) double {mustBePositive}
    
    %% Optional quantization step
    options.delta (1,1) double {mustBePositive} = Inf

    %% Optional name-value pairs for segmentation
    % Segment or not
    options.segmentation (1,1) logical = false

    % Window shape for overlap-add
    options.wtype (1,1) string = 'rect'

    % Window length for overlap-add
    options.w (1,1) double = 4096
    
    % Window shift for overlap-add
    options.a (1,1) double = 2048

    %% Optional name-value pairs for outputs
    % Save the solution during iterations
    options.saveall (1,1) logical = false

    % Plot additional linesearch graphs
    options.plotLS (1,1) logical = false
    
    % Print current iteration
    options.verbose (1,1) logical = false

    %% Optional name-value pairs for acceleration
    % Accelerate the coefficient estimation
    options.coefaccel (1,1) logical = false
    
    % Extrapolate the coefficients
    options.coefextra (1,1) logical = false
    
    % Accelerate the signal estimation
    options.sigaccel (1,1) logical = false
    
    % Extrapolate the signal
    options.sigextra (1,1) logical = false
    
    % Search for the coefficients and the signal using the linesearch strategy
    % and the current and previous solutions
    options.linesearch (1,1) logical = false

    %% Optional name-value pairs for Douglas–Rachford
    % Number of iterations for the solver of the subproblems;
    % set a number or a vector;
    % for progressive iteration, set DRmaxit=round(logspace(2, 3, maxit))'
    options.DRmaxit (1,:) double {mustBePositive} = 1000
    
    % How to build the matrices needed for the subproblems
    options.mat (1,1) string {mustBeMember(options.mat, ["toeplitz", "xcorr", "conv", "fft"])} = "toeplitz"
    
    % Use the Cholesky decomposition to substitute for the multiplication
    % with matrix inversion in the non-accelerated DR algorithm
    options.decompose (1,1) logical = true
    
    % Parameter of the DR algorithm for coefficient estimation
    options.gammaC (1,1) double {mustBePositive} = 0.1
    
    % Parameter of the DR algorithm for signal estimation
    options.gammaS (1,1) double {mustBePositive} = 10
    
end

% save the parsed results to nice variables
a            = options.a;
coefaccel    = options.coefaccel;
coefextra    = options.coefextra;
decompose    = options.decompose;
delta        = options.delta;
DRmaxit      = options.DRmaxit;
gammaC       = options.gammaC;
gammaS       = options.gammaS;
linesearch   = options.linesearch;
mat          = options.mat;
plotLS       = options.plotLS;
saveall      = options.saveall;
segmentation = options.segmentation;
sigaccel     = options.sigaccel;
sigextra     = options.sigextra;
verbose      = options.verbose;
w            = options.w;
wtype        = options.wtype;

% update plotLS
plotLS = linesearch && plotLS && ~segmentation;

% compute objective?
computeobj = nargout > 1 || linesearch;

% set the progression of the Douglas–Rachford iterations
maxits = DRmaxit(:) .* ones(maxit, 1);

% estimate delta from data if necessary
if strcmpi(method, "dequantization") && delta == Inf
    quant_levels = unique(signal);
    quant_steps = diff(quant_levels);
    delta = median(quant_steps);
    if verbose
        fprintf("Estimated quantization step: %.3f\n", delta)
    end
end

% if desired, initialize the figures for line search plotting
if plotLS
    fig1 = figure("visible", "off"); % objective function
    fig2 = figure("visible", "off"); % second derivative
    tls1 = tiledlayout(fig1, "flow");
    tls2 = tiledlayout(fig2, "flow");
else
    fig1 = [];
    fig2 = [];
    tls1 = [];
    tls2 = [];
end

%% padding the signal
if ~segmentation
    w = length(signal);
    a = length(signal);
    wtype = 'rect';
end

% L is divisible by a and minimum amount of zeros equals gl (window length)
% zeros will be appended to avoid periodization of nonzero samples
L    = ceil(length(signal)/a)*a + (ceil(w/a)-1)*a;
S    = L/a; % number of signal segments
data = [signal;  zeros(L-length(signal), 1)];
if ~strcmpi(method, "dequantization")
    masks.R = [masks.R; true(L-length(signal), 1)];
    masks.U = [masks.U; false(L-length(signal), 1)];
    masks.L = [masks.L; false(L-length(signal), 1)];
end

%% initializing the output arrays
mrestored = zeros(w, maxit, S);
objective = NaN(maxit, S);
times = NaN(maxit, S);

%% construction of analysis and synthesis windows
if strcmpi(wtype, 'rect')
    gana = ones(w, 1);
    if segmentation
        gsyn = gabwin('hann', a, w, L);
        gsyn = fftshift(gsyn);
        gsyn = normalize(gsyn, 'peak');
        gsyn = gsyn * (2 * a / w);
    else
        gsyn = gana;
    end
else
    gana = gabwin(char(wtype), a, w, L);
    gana = normalize(gana, 'peak'); % peak-normalization of the analysis window
    gana = fftshift(gana);
    gsyn = gabdual(gana, a, w)*w; % computing the synthesis window
end

% multiply the quantization step
if strcmpi(method, "dequantization")
    delta = delta * gana;
end

%% segment settings
mdata = NaN(w, S);
mR = true(w, S);
mL = false(w, S);
mU = false(w, S);
for s = 1:S
    % defining the indices of the current block
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    
    % defining the segment data and masks
    mdata(:, s) = data(indices) .* gana;
    if ~strcmpi(method, "dequantization")
        mR(:, s) = masks.R(indices);
        mL(:, s) = masks.L(indices);
        mU(:, s) = masks.U(indices);
    end
end
maxits = repmat(maxits, 1, S);

% parfor s = 1:S
for s = 1:S

    if verbose
        fprintf('Processing segment %d of %d...\n', s, S)
    end
    smasks = struct('R', mR(:, s), 'L', mL(:, s), 'U', mU(:, s));

    %% initialization
    solution = mdata(:, s);
    N = length(solution); % notational convenience, this is equal to w
    oldcoef = zeros(p+1, 1); % auxiliary variable for the extrapolation
    oldsolution = zeros(N, 1); % auxiliary variable for the extrapolation
    
    % define some useful functions
    soft = @(x, t) sign(x).*max(abs(x)-t, 0);
    if strcmpi(method, "declipping")
        proj = @(x) mdata(:, s) .* smasks.R + ...
                    max(x, mdata(:, s)) .* smasks.U + ...
                    min(x, mdata(:, s)) .* smasks.L;
    elseif strcmpi(method, "dequantization")
        proj = @(x) proj_quant(x, mdata(:, s), delta);
    else
        proj = @(x) x;
    end
    if length(lambda) > 1 && lambda(2) < Inf
        Q = @(x, c) 0.5*norm(fft(c, N+p).*fft(x, N+p))^2 / (N+p)...
            + lambda(1)*norm(c, 1)...
            + lambda(2)*0.5*norm(x-proj(x))^2;
    else
        Q = @(x, c) 0.5*norm(fft(c, N+p).*fft(x, N+p))^2 / (N+p)...
            + lambda(1)*norm(c, 1);
    end
    
    % prepare some matrices for the inpainting / glp case
    if strcmpi(method, "inpainting") || strcmpi(method, "glp")
        indmiss = find(~smasks.R);
        indobs  = find(smasks.R);
        IAA     = abs(repmat(indmiss, 1, N)-repmat(1:N, length(indmiss), 1));
        IAA1    = IAA <= p;
    else
        indmiss = [];
        indobs = [];
        IAA = [];
        IAA1 = [];
    end
    
    %% main iteration
    str = "";
    t = tic;
    for i = 1:maxit
    
        if verbose
            fprintf(repmat('\b', 1, strlength(str)))
            str = sprintf("iteration %d of %d", i, maxit);
            fprintf(str)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           AR model estimation                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if lambda(1) == 0        
            coef = lpc(solution, p)'; 
        else
            % prepare the matrix XX
            if strcmpi(mat, "toeplitz")
                X  = toeplitz([solution; zeros(p, 1)], [solution(1), zeros(1, p)]);
                XX = X'*X;
            elseif strcmpi(mat, "xcorr")
                x  = xcorr(solution, p)';
                XX = spdiags(ones(p+1, 1)*x, -p:p, p+1, p+1);
            elseif strcmpi(mat, "conv")
                x  = conv(solution, flip(solution))';
                XX = spdiags(ones(p+1, 1)*x(N-p:N+p), -p:p, p+1, p+1);
            elseif strcmpi(mat, "fft")
                x  = ifft(fft([solution' zeros(1, p)]) .* fft([flip(solution)' zeros(1, p)]));
                XX = spdiags(ones(p+1, 1)*x(N-p:N+p), -p:p, p+1, p+1);
            else
                XX = [];
            end
            
            % check the positive definiteness using chol
            try
                dXX = decomposition(eye(p+1) + gammaC*XX, "chol");
            catch
                if verbose
                    warning("Stopped on signal matrix decomposition in iteration %d.", i)
                end
                break
            end
            
            % check the conditioning
            if isIllConditioned(dXX)
                if verbose
                    warning("Stopped on signal matrix ill conditioning in iteration %d.", i)
                end
                break
            end
            
            if coefaccel
                % set the parameters of the accelerated Douglas–Rachford
                % algorithm
                P     = @(a) lambda(1)*norm(a, 1);
                proxP = @(a, t) [1; soft(a(2:end), lambda(1)*t)]; % proximal operator of the penalty
    
                % solve the task
                coef = DouglasRachfordA(zeros(N+p, 1), solution, P, proxP, p+1, maxits(i, s), "alpha", gammaC);
            else                
                % set the parameters of the model
                f  = @(a) 0.5*norm(fft(solution, N+p).*fft(a, N+p))^2 / (N+p);
                if decompose
                    prox_f = @(a, t) dXX\a;
                else
                    invmat = inv(eye(p+1) + gammaC*XX);
                    prox_f = @(a, t) invmat*a; %#ok<MINV>
                end
                g      = @(a) lambda(1)*norm(a(2:end), 1);
                prox_g = @(a, t) [1; soft(a(2:end), lambda(1)*t)]; 
    
                % solve the task
                coef = DouglasRachford(f, g, prox_f, prox_g, p+1, ...
                    "lambda", 1, "gamma", gammaC, "maxit", maxits(i, s), ...
                    "startpoint", [1; zeros(p, 1)], "tol", -Inf);
            end
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          AR model extrapolation                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if coefextra && i > 1
            % this formula is heuristic
            coef = coef + 2*(maxit-i)/(maxit-1) * (coef - oldcoef);
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            signal estimation                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmpi(method, "inpainting") || strcmpi(method, "glp")
            AA       = zeros(size(IAA));
            b        = coef'*hankel(coef', [coef(end), zeros(1, p)]);
            AA(IAA1) = b(IAA(IAA1)+1);
            [R, er]  = chol(AA(:, indmiss));
            if er
                if verbose
                    warning("Stopped on AA matrix decomposition in iteration %d.", i)
                end
                break
            else
                solution(~smasks.R) = -R\(R'\(AA(:, indobs)*mdata(indobs, s)));
            end
        else
            % prepare the matrix AA
            if strcmpi(mat, "toeplitz")
                A  = toeplitz([coef; zeros(N-1, 1)], [coef(1), zeros(1, N-1)]);
                AA = A'*A;
            elseif strcmpi(mat, "xcorr")
                b  = xcorr(coef, p)';
                AA = spdiags(ones(N, 1)*b, -p:p, N, N);
            elseif strcmpi(mat, "conv")
                b  = conv(coef, flip(coef))';
                AA = spdiags(ones(N, 1)*b, -p:p, N, N);
            elseif strcmpi(mat, "fft")
                b  = ifft(fft([coef' zeros(1, p)]) .* fft([flip(coef') zeros(1, p)]));
                AA = spdiags(ones(N, 1)*b, -p:p, N, N);
            else
                AA = [];
            end
    
            % check the positive definiteness using chol
            try
                dAA = decomposition(eye(N) + gammaS*AA, "chol");
            catch
                if verbose
                    warning("Stopped on coef matrix decomposition in iteration %d.", i)
                end
                break
            end
            
            % check the conditioning
            if isIllConditioned(dAA)
                if verbose
                    warning("Stopped on coef matrix ill conditioning in iteration %d.", i)
                end
                break
            end
    
            if sigaccel
                % set the parameters of the accelerated Douglas–Rachford algorithm
                if length(lambda) > 1  && lambda(2) < Inf
                    P     = @(x) lambda(2)*0.5*norm(x-proj(x))^2;
                    proxP = @(x, t) lambda(2)*t/(lambda(2)*t+1)*proj(x) + 1/(lambda(2)*t+1)*x;
                else
                    P     = @(x) 0; % this is actually the indicator function
                    proxP = @(x, t) proj(x);
                end
    
                % solve the task
                solution = DouglasRachfordA(zeros(N+p, 1), coef, P, proxP, N, maxits(i, s), "alpha", gammaS);
            else   
                % set the parameters of the model
                f = @(x) 0.5*norm(fft(coef, N+p).*fft(x, N+p))^2 / (N+p);
                if decompose
                    prox_f = @(x, t) dAA\x;
                else
                    invmat = inv(eye(N) + gammaS*AA);
                    prox_f = @(x, t) invmat*x; %#ok<MINV>
                end
                if length(lambda) > 1  && lambda(2) < Inf
                    g      = @(x) lambda(2)*0.5*norm(x-proj(x))^2;
                    prox_g = @(x, t) lambda(2)*t/(lambda(2)*t+1)*proj(x) + 1/(lambda(2)*t+1)*x;
                else
                    g      = @(x) 0;
                    prox_g = @(x, t) proj(x);
                end
    
                % solve the task
                solution = DouglasRachford(f, g, prox_f, prox_g, N, ...
                    "lambda", 1, "gamma", gammaS, "maxit", maxits(i, s), ...
                    "startpoint", solution, "tol", -Inf);
            end
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    signal rectification (for GLP only)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmpi(method, "glp") && i < maxit        
            % signal segments to be checked
            solutionU = solution(smasks.U);
            solutionL = solution(smasks.L);
            signalU   = mdata(smasks.U, s);
            signalL   = mdata(smasks.L, s);
            
            % indicators of the samples to be rectified
            IU = solutionU < signalU;
            IL = solutionL > signalL;
                    
            % rectification
            solutionU(IU)      = signalU(IU) + abs(solutionU(IU) - signalU(IU));
            solutionL(IL)      = signalL(IL) - abs(solutionL(IL) - signalL(IL));
            solution(smasks.U) = solutionU;
            solution(smasks.L) = solutionL;
        end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           signal extrapolation                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if sigextra && i > 1
            % this formula is heuristic
            solution = solution + (maxit-i)/(maxit-1) * (solution - oldsolution);
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                linesearch                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if linesearch && i > 1 && i < maxit
            
            % search for the solution and coefs on the line given the previous
            % and current solution and previous and current coefs
            steps = logspace(-4, 2, 100);
            QQ = zeros(length(steps), 1);
            for j = 1:length(steps)
                QQ(j) = Q(solution + steps(j)*(solution-oldsolution), ...
                    coef + steps(j)*(coef-oldcoef));
            end
            
            % plot the objective function and the second derivative along the
            % line
            if plotLS
                x1 = solution;
                x2 = x1 - oldsolution;
                A1 = toeplitz([coef; zeros(N-1, 1)], [coef(1), zeros(1, N-1)]);
                A2 = A1 - toeplitz([oldcoef; zeros(N-1, 1)], [oldcoef(1), zeros(1, N-1)]);
                C0 = x1'*(A1'*A1)*x1;
                C1 = 2*(x2'*A1' + x1'*A2')*A1*x1;
                C2 = x1'*(A2'*A2)*x1 + x2'*(A1'*A1)*x2 + 2*x2'*(A2'*A1 + A1'*A2)*x1;
                C3 = 2*(x1'*A2' + x2'*A1')*A2*x2;
                C4 = x2'*(A2'*A2)*x2;
                pol1 = @(t) 0.5*(C4*t.^4 + C3*t.^3 + C2*t.^2 + C1*t + C0);
                pol2 = @(t) 6*C4*t.^2 + 3*C3*t + C2;
    
                % plot the objective function
                nexttile(tls1)
                h1 = loglog(steps, pol1(steps));
                hold on
                h2 = loglog(steps, QQ);
                [m1, I] = min(QQ);
                x1 = steps(I);
                Q2min = @(step) Q(solution + step*(solution-oldsolution), ...
                        coef + step*(coef-oldcoef)); 
                [x2, m2] = fminbnd(Q2min, steps(1), steps(end));
                m0 = Q(solution, coef); % objective with no extrapolation
                title(sprintf(...
                    "iteration %d\n" + ...
                    "min. point from samples: [%.3e, %.3e]\n" + ...
                    "using fminbnd: [%.3e, %.3e]", ...
                    i, x1, m1-m0, x2, m2-m0))
                xline(m1, ":k")
                xlabel("step")
                ylabel("objective")
                legend([h1, h2], "AR part", "whole objective")
    
                % plot the secont derivative
                nexttile(tls2)
                numerical = gradient(gradient(QQ));
                polynomial = pol2(steps);
                h1 = loglog(steps(polynomial > 0), polynomial(polynomial > 0));
                hold on
                h2 = loglog(steps(numerical > 0), numerical(numerical > 0));
                title(sprintf("iteration %d", i))
                xline(m1, ":k")
                xlabel("step")
                ylabel("2nd derivative (if positive)")
                xlim(steps([1, end]))
                legend([h1, h2], "polynomial (only AR part)", "numerical (whole objective)")
            end
    
            % update the solution
            [minimum, index] = min(QQ);
            if minimum < Q(solution, coef)
                solution = solution + steps(index)*(solution-oldsolution);
                coef = coef + steps(index)*(coef-oldcoef);
            end
            
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       solution and objective update                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % update the solution
        mrestored(:, i, s) = solution;
    
	    % compute the objective value
        if computeobj
            objective(i, s) = Q(solution, coef);
        end
        
        % save the current solution for the next iteration
        oldcoef = coef;
        oldsolution = solution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                time update                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        times(i, s) = toc(t);
    
    end
    
    if verbose
        fprintf(repmat('\b', 1, strlength(str)))
    end

end

%% overlap-add
restored = zeros(L, maxit);
for s = 1:S
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    restored(indices, :) = restored(indices, :) + mrestored(:, :, s).*repmat(gsyn, 1, maxit);
end

%% cropping the solution to the original length
restored = restored(1:length(signal), :);
if ~saveall
    restored = restored(:, end);
    objective = sum(objective, 2);
    times = sum(times, 2);
end

%% rename the figures
if plotLS
    title(tls1, "lineaserach in " + method + ", " + ...
        "showing also optimal step and objective change " + ...
        "(negative means improvement)")
    title(tls2, "linesearch in " + method + ", derivatives")

    set(fig1, "WindowState", "maximized")
    set(fig1, "visible", "on")
    set(fig2, "WindowState", "maximized")
    set(fig2, "visible", "on")
end

end