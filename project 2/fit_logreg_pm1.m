function [beta, step] = fit_logreg_pm1(Phi, b, beta0)
%FIT_LOGREG_PM1 fit logistische regressie met labels in {-1,1} via GD

if nargin < 3 || isempty(beta0)
    beta0 = zeros(size(Phi, 2), 1);
end

f = @(x) logreg_cost_pm1(x, Phi, b);
df = @(x) logreg_grad_pm1(x, Phi, b);

% Param: step_size, tol, max_iter, armijo_c, backtrack_rho
[beta, step] = GD(f, df, beta0, 1, 1e-6, 15000, 1e-4, 0.5);

end
