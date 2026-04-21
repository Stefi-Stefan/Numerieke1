function [yhat, nWrong] = predict_pm1(beta, Phi, ytrue)
%PREDICT_PM1 Voorspelt labels in {-1,1} op basis van logit-teken

scores = Phi * beta;
yhat = sign(scores);
yhat(yhat == 0) = 1;

if nargin >= 3 && ~isempty(ytrue)
    nWrong = sum(yhat ~= ytrue);
else
    nWrong = [];
end

end
