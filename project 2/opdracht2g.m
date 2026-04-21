% OPDRACHT 2G
% Leave-One-Out Cross-Validation voor n = 1..12

clearvars -except quickMode; clc; close all;

% quickMode=true voor snellere test (minder n-waarden)
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
if quickMode
    nVals = 1:3;
else
    nVals = 1:12;
end

% LOOCV gebruikt veel model-fits; gebruik iets lossere GD-instellingen
% dan in 2(d) om de rekentijd beheersbaar te houden
if quickMode
    gdTol = 1e-5;
    gdMaxIter = 1500;
else
    gdTol = 1e-5;
    gdMaxIter = 4000;
end

cvLoo = zeros(size(nVals));
avgSteps = zeros(size(nVals));

for ii = 1:numel(nVals)
    n = nVals(ii);

    % precompute volledige featurematrix voor deze n
    PhiAll = build_poly_features(x1, x2, n);

    nWrong = 0;
    stepSum = 0;
    betaInit = zeros(2*n + 1, 1);

    for k = 1:N
        mask = true(N, 1);
        mask(k) = false;

        PhiTr = PhiAll(mask, :);
        bTr = b(mask);
        PhiTe = PhiAll(k, :);

        f = @(x) logreg_cost_pm1(x, PhiTr, bTr);
        df = @(x) logreg_grad_pm1(x, PhiTr, bTr);
        [beta, step] = GD(f, df, betaInit, 1, gdTol, gdMaxIter, 1e-4, 0.5);
        betaInit = beta; % warm start voor volgende fold
        stepSum = stepSum + step;

        yhat = predict_pm1(beta, PhiTe);
        if yhat ~= b(k)
            nWrong = nWrong + 1;
        end
    end

    cvLoo(ii) = nWrong / N;
    avgSteps(ii) = stepSum / N;

    fprintf('n = %2d: CV_LOO = %.4f (gem. GD stappen = %.1f)\n', n, cvLoo(ii), avgSteps(ii));
end

[bestErr, idxBest] = min(cvLoo);
bestN = nVals(idxBest);

fprintf('Beste model volgens LOOCV: n = %d met CV_LOO = %.4f\n', bestN, bestErr);

if ~quickMode
    outDir = fullfile(pwd, 'Figuren');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fig = figure('Color', 'w');
    plot(nVals, cvLoo, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
    hold on;
    plot(bestN, bestErr, 'rs', 'MarkerFaceColor', 'r');
    hold off;

    xlabel('Modelgraad n');
    ylabel('LOOCV-kruisvalidatiefout (ratio)');
    title('Opgave 2(g): LOOCV-fout in functie van n');
    legend({'CV^{LOO}(n)', sprintf('minimum: n=%d', bestN)}, 'Location', 'best');
    grid on;

    saveas(fig, fullfile(outDir, 'opgave2g_loocv.eps'), 'epsc');
    saveas(fig, fullfile(outDir, 'opgave2g_loocv.png'));
end
