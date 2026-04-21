% OPDRACHT 2I
% Herhaal 10-fold cross-validatie vijf keer na het schudden van de data

clearvars -except quickMode; clc; close all;

if ~exist('quickMode', 'var')
    quickMode = false;
end

S = load('DatasetCV.mat');
x1 = S.x(:);
x2 = S.y(:);
b  = S.cat(:);

if ~all(ismember(unique(b), [-1, 1]))
    error('Labels in DatasetCV.mat moeten in {-1,1} liggen.');
end

N = numel(b);
K = 10;

if quickMode
    nVals = 1:4;
else
    nVals = 1:20;
end

% lossere GD-instellingen dan in 2(d), omdat hier veel fits nodig zijn
if quickMode
    gdTol = 1e-5;
    gdMaxIter = 2000;
else
    gdTol = 1e-5;
    gdMaxIter = 4500;
end

nRuns = 5;
cvAll = zeros(nRuns, numel(nVals));
bestN = zeros(nRuns, 1);
bestErr = zeros(nRuns, 1);

for r = 1:nRuns
    % randomise/schud de originele data per run
    rng(100 + r, 'twister');
    perm = randperm(N);
    x1s = x1(perm);
    x2s = x2(perm);
    bs  = b(perm);

    foldId = mod(0:N-1, K) + 1;

    for ii = 1:numel(nVals)
        n = nVals(ii);
        PhiAll = build_poly_features(x1s, x2s, n);

        wrong = 0;

        for k = 1:K
            idxTe = (foldId == k);
            idxTr = ~idxTe;

            PhiTr = PhiAll(idxTr, :);
            PhiTe = PhiAll(idxTe, :);
            bTr = bs(idxTr);
            bTe = bs(idxTe);

            beta0 = zeros(size(PhiTr, 2), 1);
            f = @(x) logreg_cost_pm1(x, PhiTr, bTr);
            df = @(x) logreg_grad_pm1(x, PhiTr, bTr);
            beta = GD(f, df, beta0, 1, gdTol, gdMaxIter, 1e-4, 0.5);

            yhat = predict_pm1(beta, PhiTe);
            wrong = wrong + sum(yhat ~= bTe);
        end

        cvAll(r, ii) = wrong / N;
    end

    [bestErr(r), idxBest] = min(cvAll(r, :));
    bestN(r) = nVals(idxBest);

    fprintf('Run %d: beste n = %d, CV = %.4f\n', r, bestN(r), bestErr(r));
end

outDir = fullfile(pwd, 'Figuren');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fig = figure('Color', 'w');
hold on;
cols = lines(nRuns);
for r = 1:nRuns
    plot(nVals, cvAll(r, :), '-', 'Color', cols(r, :), 'LineWidth', 1.3);
end
hold off;

xlabel('Modelgraad n');
ylabel('10-fold kruisvalidatiefout (ratio)');
title('Opgave 2(i): 10-fold CV na shuffelen van de data');
grid on;

legTxt = arrayfun(@(r) sprintf('shuffle %d', r), 1:nRuns, 'UniformOutput', false);
legend(legTxt, 'Location', 'best');

saveas(fig, fullfile(outDir, 'opgave2i_kfcv10_herhalingen.eps'), 'epsc');
saveas(fig, fullfile(outDir, 'opgave2i_kfcv10_herhalingen.png'));
