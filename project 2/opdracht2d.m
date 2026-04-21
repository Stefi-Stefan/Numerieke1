% OPDRACHT 2D
% Modelleert de dataset met logistische regressie voor n = 0..5
% voor elk n wordt een figuur gemaakt met datapunten en p=1/2-scheiding

clear; clc; close all;

S = load('DatasetCV.mat');
x1 = S.x(:);
x2 = S.y(:);
b  = S.cat(:);

if ~all(ismember(unique(b), [-1, 1]))
    error('Labels in DatasetCV.mat moeten in {-1,1} liggen.');
end

outDir = fullfile(pwd, 'Figuren');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% raster voor decision boundary
nGrid = 220;
x1g = linspace(min(x1), max(x1), nGrid);
x2g = linspace(min(x2), max(x2), nGrid);
[X1, X2] = meshgrid(x1g, x2g);

for n = 0:5
    Phi = build_poly_features(x1, x2, n);
    [beta, step] = fit_logreg_pm1(Phi, b);

    [~, nWrong] = predict_pm1(beta, Phi, b);
    fprintf('n = %d: fout geclassificeerd = %d (GD stappen = %d)\n', n, nWrong, step);

    % evaluatie op raster: p=1/2 komt overeen met logit = 0.
    PhiGrid = build_poly_features(X1(:), X2(:), n);
    logitGrid = reshape(PhiGrid * beta, size(X1));

    fig = figure('Color', 'w');
    hold on;

    idxPos = (b == 1);
    idxNeg = (b == -1);

    scatter(x1(idxPos), x2(idxPos), 26, [0.05 0.35 0.90], 'filled');
    scatter(x1(idxNeg), x2(idxNeg), 26, [0.85 0.15 0.10], 'filled');

    hasBoundary = (min(logitGrid(:)) < 0) && (max(logitGrid(:)) > 0);
    if hasBoundary
        contour(X1, X2, logitGrid, [0 0], 'k', 'LineWidth', 1.6);
    else
        text(mean(x1g), max(x2g), 'geen p=1/2-scheiding in plotgebied', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    xlabel('x_1');
    ylabel('x_2');
    title(sprintf('Logistische regressie, n = %d, fout = %d', n, nWrong));
    if hasBoundary
        legend({'klasse 1', 'klasse -1', 'p = 1/2'}, 'Location', 'best');
    else
        legend({'klasse 1', 'klasse -1'}, 'Location', 'best');
    end
    grid on;
    axis tight;
    hold off;

    % EPS voor verslag, PNG voor preview in de file explorer
    saveas(fig, fullfile(outDir, sprintf('opgave2d_n%d.eps', n)), 'epsc');
    saveas(fig, fullfile(outDir, sprintf('opgave2d_n%d.png', n)));
end
