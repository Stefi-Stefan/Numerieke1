function J = logreg_cost_pm1(beta, Phi, b)
%LOGREG_COST_PM1 gemiddelde kost voor b in {-1,1}

z = Phi * beta;
t = -b .* z;

%numeriek stabiele softplus: log(1 + exp(t)).
J = mean(max(t, 0) + log1p(exp(-abs(t))));

end
