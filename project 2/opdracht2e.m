% OPDRACHT 2E
% Random split D = Tr U Te (even groot) en CV-fout voor n = 1..20

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

% reproduceerbare random split
rng(1, 'twister');
perm = randperm(N);
nTrain = N/2;
idxTr = perm(1:nTrain);
idxTe = perm(nTrain+1:end);

x1Tr = x1(idxTr); x2Tr = x2(idxTr); bTr = b(idxTr);
x1Te = x1(idxTe); x2Te = x2(idxTe); bTe = b(idxTe);

nVals = 1:20;
cvErr = zeros(size(nVals));
gdSteps = zeros(size(nVals));

for i = 1:numel(nVals)
    n = nVals(i);

    PhiTr = build_poly_features(x1Tr, x2Tr, n);
    PhiTe = build_poly_features(x1Te, x2Te, n);

    [beta, step] = fit_logreg_pm1(PhiTr, bTr);
    gdSteps(i) = step;

    [~, nWrong] = predict_pm1(beta, PhiTe, bTe);
    cvErr(i) = nWrong / numel(bTe);
end

[bestErr, bestIdx] = min(cvErr);
bestN = nVals(bestIdx);

for i = 1:numel(nVals)
    fprintf('n = %2d: CV = %.4f (GD stappen = %d)\n', nVals(i), cvErr(i), gdSteps(i));
end
fprintf('Beste model in deze split: n = %d met CV = %.4f\n', bestN, bestErr);

outDir = fullfile(pwd, 'Figuren');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fig = figure('Color', 'w');
plot(nVals, cvErr, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
hold on;
plot(bestN, bestErr, 'rs', 'MarkerFaceColor', 'r');
hold off;

xlabel('Modelgraad n');
ylabel('Kruisvalidatiefout (ratio)');
title('Opgave 2(e): CV-fout voor random train/test split');
grid on;
legend({'CV(n)', sprintf('minimum: n=%d', bestN)}, 'Location', 'best');

saveas(fig, fullfile(outDir, 'opgave2e_cv_random.eps'), 'epsc');
saveas(fig, fullfile(outDir, 'opgave2e_cv_random.png'));
