function [x, step] = GD(func, Dfunc, x0, step_size, tol, max_iter, armijo_c, backtrack_rho)
%GD  Gradient Descent met backtracking line search
%   [x, step] = GD(func, Dfunc, x0, ...)
%   - func:  function handle naar f(x)
%   - Dfunc: function handle naar grad f(x)
%   - x0:    initiële gok
%
%   Optionele parameters (met defaults):
%   step_size     initiële stapgrootte (default 1)
%   tol           tolerantie voor stopcriteria (default 1e-6)
%   max_iter      maximum aantal iteraties (default 5000)
%   armijo_c      Armijo-parameter c in (0,1) (default 1e-4)
%   backtrack_rho reductiefactor rho in (0,1) (default 0.5)

if nargin == 3
    step_size = 1;
    tol = 1e-6;
    max_iter = 5000;
    armijo_c = 1e-4;
    backtrack_rho = 0.5;
elseif nargin < 8
    if nargin < 4 || isempty(step_size), step_size = 1; end
    if nargin < 5 || isempty(tol), tol = 1e-6; end
    if nargin < 6 || isempty(max_iter), max_iter = 5000; end
    if nargin < 7 || isempty(armijo_c), armijo_c = 1e-4; end
    if nargin < 8 || isempty(backtrack_rho), backtrack_rho = 0.5; end
end

x = x0;
step = 0;

for k = 1:max_iter
    g = Dfunc(x);

    % Stop bij kleine gradiëntnorm.
    if norm(g, 2) <= tol
        step = k - 1;
        return;
    end

    alpha = step_size;
    fx = func(x);
    g2 = g(:).' * g(:);

    % Backtracking op basis van Armijo-voorwaarde.
    while func(x - alpha * g) > fx - armijo_c * alpha * g2
        alpha = backtrack_rho * alpha;

        % Veiligheidsstop bij extreem kleine stap (machine precisie)
        if alpha < 1e-16
            step = k - 1;
            return;
        end
    end

    x_new = x - alpha * g;

    % Stop bij kleine relatieve iteratieverandering
    if norm(x_new - x, 2) <= tol * (1 + norm(x, 2))
        x = x_new;
        step = k;
        return;
    end

    x = x_new;
    step = k;
end
end
