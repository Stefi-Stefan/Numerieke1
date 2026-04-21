% OPDRACHT 2H
% 10-fold cross-validatie voor n = 1..20

clearvars -except quickMode; clc; close all;

% quickMode=true voor snellere test met minder n-waarden
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

% reproduceerbare random verdeling van de data over 10 vouwen
rng(10, 'twister');
perm = randperm(N);
foldId = mod(0:N-1, K) + 1;
foldId = foldId(perm);

cvK = zeros(size(nVals));
avgSteps = zeros(size(nVals));

for ii = 1:numel(nVals)
    n = nVals(ii);
    PhiAll = build_poly_features(x1, x2, n);

    wrong = 0;
    stepSum = 0;

    for k = 1:K
        idxTe = (foldId == k);
        idxTr = ~idxTe;

        PhiTr = PhiAll(idxTr, :);
        PhiTe = PhiAll(idxTe, :);
        bTr = b(idxTr);
        bTe = b(idxTe);

        beta0 = zeros(size(PhiTr, 2), 1);
        f = @(x) logreg_cost_pm1(x, PhiTr, bTr);
        df = @(x) logreg_grad_pm1(x, PhiTr, bTr);
        [beta, step] = GD(f, df, beta0, 1, gdTol, gdMaxIter, 1e-4, 0.5);

        stepSum = stepSum + step;

        yhat = predict_pm1(beta, PhiTe);
        wrong = wrong + sum(yhat ~= bTe);
    end

    cvK(ii) = wrong / N;
    avgSteps(ii) = stepSum / K;

    fprintf('n = %2d: CV_10-fold = %.4f (gem. GD stappen = %.1f)\n', n, cvK(ii), avgSteps(ii));
end

[bestErr, idxBest] = min(cvK);
bestN = nVals(idxBest);

fprintf('Beste model volgens 10-fold CV: n = %d met CV = %.4f\n', bestN, bestErr);

if ~quickMode
    outDir = fullfile(pwd, 'Figuren');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fig = figure('Color', 'w');
    plot(nVals, cvK, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
    hold on;
    plot(bestN, bestErr, 'rs', 'MarkerFaceColor', 'r');
    hold off;

    xlabel('Modelgraad n');
    ylabel('10-fold kruisvalidatiefout (ratio)');
    title('Opgave 2(h): 10-fold CV-fout in functie van n');
    legend({'CV^{10-fold}(n)', sprintf('minimum: n=%d', bestN)}, 'Location', 'best');
    grid on;

    saveas(fig, fullfile(outDir, 'opgave2h_kfcv10.eps'), 'epsc');
    saveas(fig, fullfile(outDir, 'opgave2h_kfcv10.png'));
end
