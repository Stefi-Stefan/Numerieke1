% OPDRACHT 2F
% Herhaalt opgave 2(e) vijf keer met nieuwe random splits
% Deze figuur bevat dus 6 CV-curves

clear; clc; close all;

S = load('DatasetCV.mat');
x1 = S.x(:);
x2 = S.y(:);
b  = S.cat(:);

if ~all(ismember(unique(b), [-1, 1]))
    error('Labels in DatasetCV.mat moeten in {-1,1} liggen.');
end

N = numel(b);
if mod(N, 2) ~= 0
    error('Voor deze opgave wordt een even aantal datapunten verwacht.');
end

nVals = 1:20;
nRuns = 6; % 1 run uit 2(e) + 5 extra runs
cvAll = zeros(nRuns, numel(nVals));
bestN = zeros(nRuns, 1);
bestErr = zeros(nRuns, 1);

for r = 1:nRuns
    rng(r, 'twister');
    perm = randperm(N);

    nTrain = N/2;
    idxTr = perm(1:nTrain);
    idxTe = perm(nTrain+1:end);

    x1Tr = x1(idxTr); x2Tr = x2(idxTr); bTr = b(idxTr);
    x1Te = x1(idxTe); x2Te = x2(idxTe); bTe = b(idxTe);

    for i = 1:numel(nVals)
        n = nVals(i);
        PhiTr = build_poly_features(x1Tr, x2Tr, n);
        PhiTe = build_poly_features(x1Te, x2Te, n);

        beta = fit_logreg_pm1(PhiTr, bTr);
        [~, nWrong] = predict_pm1(beta, PhiTe, bTe);
        cvAll(r, i) = nWrong / numel(bTe);
    end

    [bestErr(r), idxMin] = min(cvAll(r, :));
    bestN(r) = nVals(idxMin);

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

xlabel('Modelgraad n');
ylabel('Kruisvalidatiefout (ratio)');
title('Opgave 2(f): CV-fout voor 6 random train/test splits');
grid on;

legTxt = arrayfun(@(r) sprintf('split %d', r), 1:nRuns, 'UniformOutput', false);
legend(legTxt, 'Location', 'best');

hold off;

saveas(fig, fullfile(outDir, 'opgave2f_cv_random_herhalingen.eps'), 'epsc');
saveas(fig, fullfile(outDir, 'opgave2f_cv_random_herhalingen.png'));
