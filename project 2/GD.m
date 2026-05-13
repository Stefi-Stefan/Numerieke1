function [x, step] = GD(func, Dfunc, x0, step_size, eps, max_iter, gamma, q)
%GD  Gradient Descent met backtracking line search
%   [x, step] = GD(func, Dfunc, x0, ...)
%   - x: de oplossing
%   - step: de stap waarin het algoritme stopt
%
%   - func:  function handle naar f(x)
%   - Dfunc: function handle naar grad f(x)
%   - x0:    initiële gok
%
%   Optionele parameters (met defaults):
%   step_size     initiële stapgrootte (default 1)
%   eps           tolerantie voor stopcriteria (default 1e-6)
%   max_iter      maximum aantal iteraties (default 5000)
%   gamma      Armijo-parameter gamma in (0,1) (default 1e-4)
%   q reductiefactor in (0,1) (default 0.5)

if nargin == 3
    step_size = 1;
    eps = 1e-6;
    max_iter = 5000;
    gamma = 1e-4;
    q = 0.5;
else
    if nargin < 4 || isempty(step_size), step_size = 1; end
    if nargin < 5 || isempty(eps), eps = 1e-6; end
    if nargin < 6 || isempty(max_iter), max_iter = 5000; end
    if nargin < 7 || isempty(gamma), gamma = 1e-4; end
    if nargin < 8 || isempty(q), q = 0.5; end
end

x = x0;
step = 0;

for k = 1:max_iter
    g = Dfunc(x);

    % stoppen als 2-norm vd gradient kleiner is dan tolerantie,
    % dan is gradient heel klein en ongeveer lokaal minimum bereikt
    if norm(g, 2) <= eps
        step = k - 1;
        return;
    end

    alpha = step_size;
    fx = func(x);
    g2 = g(:).' * g(:);

    %backtracking op basis van Armijo-voorwaarde:
    % f(x^(k-1) + alpha * p^(k)) <= f(x^(k-1)) + gamma * alpha * grad(f(x^(k-1)))^T * p^(k)
    % hier is richting p^(k) = -g , dus
    % f(x - alpha * g) <= f(x) + gamma * alpha * g^T * (-g) 
    % f(x - alpha * g) <= f(x) - gamma * alpha * ||g||^2_2
    while func(x - alpha * g) > fx - gamma * alpha * g2
        alpha = q * alpha;

        % bij step met grootte machine precisie stoppen
        if alpha < 1e-16
            step = k - 1;
            return;
        end
    end

    x_new = x - alpha * g;
    % stopcriterium, bij kleine relatieve verandering van het iteratiepunt stoppen
    if norm(x_new - x, 2) <= eps * (1 + norm(x, 2))
        x = x_new;
        step = k;
        return;
    end
    x = x_new;
    step = k;
end
end
