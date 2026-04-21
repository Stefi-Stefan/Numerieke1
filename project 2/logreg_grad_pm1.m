function g = logreg_grad_pm1(beta, Phi, b)
%LOGREG_GRAD_PM1 gradiënt van logreg_cost_pm1

z = Phi * beta;
den = 1 + exp(b .* z);
g = -(Phi' * (b ./ den)) / numel(b);

end
