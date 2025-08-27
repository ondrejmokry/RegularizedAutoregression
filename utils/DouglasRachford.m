function [x_hat, obj_val, rel_norm] = DouglasRachford(f, g, prox_f, prox_g, dim, options)
% Douglas-Rachford algorithm for solving
%
%      x_hat = arg min f(x) + g(x)
%
% algorithm taken from
%     Combettes, Patrick & Pesquet, Jean-Christophe. (2009). Proximal
%     Splitting Methods in Signal Processing. Fixed-Point Algorithms for
%     Inverse Problems in Science and Engineering. 49.
%     10.1007/978-1-4419-9569-8_10.
%
% outputs:
%     x_hat        solution
%     obj_val      value of f(x) + g(x) during the iterations
%     rel_norm     relative norm of x wrt previous iteration
%
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

arguments
    f function_handle % function handle for f(x)
    g function_handle % function handle for g(x)
    prox_f function_handle % proximal operator of f(x), used as prox_f(arg, parameter)
    prox_g function_handle % proximal operator of g(x), used as prox_g(arg, parameter)
    dim (1,1) double {mustBePositive, mustBeInteger} % length of x
    options.lambda (1,1) double = 1 % step size
    options.gamma (1,1) double = 1 % parameter of the proximal operators
    options.startpoint (:,1) double = zeros(dim, 1) % starting point
    options.maxit (1,1) double {mustBePositive, mustBeInteger} = 500 % maximal number of iterations
    options.tol (1,1) double = 5e-4 % tolerance of the relative norm
end

% Check that the startpoint has correct length
assert(length(options.startpoint) == dim, 'Starting point y0 must have length equal to dim');

%% initialization
obj_val = Inf(options.maxit, 1);
rel_norm = Inf(options.maxit, 1);
x_old = options.startpoint;
y = options.startpoint;

%% algorithm
for i = 1:options.maxit
    x_new = prox_g(y, options.gamma);
    y = y + options.lambda*(prox_f(2*x_new - y, options.gamma) - x_new);
    if nargout > 1
        obj_val(i) = f(x_new) + g(x_new);
    end
    rel_norm(i) = norm(x_new-x_old)/norm(x_old);
    x_old = x_new;
    if  rel_norm(i) <= options.tol
        break
    end
end

%% output
x_hat = x_new;
if nargout > 1
    obj_val = obj_val(1:i-1);
    rel_norm = rel_norm(1:i-1);
end