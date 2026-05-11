% OPDRACHT 2J
% Bias en variantie van LOOCV-schatter vs. dataset-grootte
% Voor N = 10 + 20i, i = 1..19, worden 50 samples getrokken en LOOCV berekend

clear; close all;

% Model degree - gemakkelijk aan te passen
n = 3;

% Laad grote dataset
load('BigDatasetCV.mat', 'x', 'y', 'cat');
x1_big = x(:);
x2_big = y(:);
b_big = cat(:);
N_big = length(x1_big);
X_big = [x1_big, x2_big];

% Grootteparameters: N = 10 + 20i voor i = 1..19
N_values = 10 + 20 * (1:19);
num_N = length(N_values);
num_samples = 50;

mean_cv = zeros(1, num_N);
var_cv = zeros(1, num_N);

tic;
fprintf('Opdracht 2(j): %d N-waarden, %d samples per N, n=%d\n', num_N, num_samples, n);

for i_N = 1:num_N
    N = N_values(i_N);
    fprintf('N = %3d (%2d/%2d): ', N, i_N, num_N);
    
    cv_values = zeros(1, num_samples);
    
    % Genereer 50 willekeurige samples van grootte N (zonder teruglegging)
    for i_sample = 1:num_samples
        idx = randperm(N_big, N);
        X_sub = X_big(idx, :);
        b_sub = b_big(idx);
        
        % Bereken LOOCV voor deze sample
        cv_loo = compute_loocv_subset(X_sub, b_sub, n);
        cv_values(i_sample) = cv_loo;
    end
    
    mean_cv(i_N) = mean(cv_values);
    var_cv(i_N) = var(cv_values);
    
    fprintf('mean=%.4f, var=%.6f\n', mean_cv(i_N), var_cv(i_N));
end

t = toc;
fprintf('Totale tijd: %.1f seconden\n\n', t);

% Figuur 1: Gemiddelde LOOCV-fout vs. N
figure('Position', [100 100 900 600]);
plot(N_values, mean_cv, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Dataset-grootte N', 'FontSize', 12);
ylabel('Gemiddelde CV_{LOO}', 'FontSize', 12);
title(sprintf('Gemiddelde LOOCV-fout vs. dataset-grootte (n=%d)', n), 'FontSize', 13);
grid on;
xlim([min(N_values)-10, max(N_values)+10]);

saveas(gcf, 'Figuren/opgave2j_mean_cv.eps', 'epsc');
saveas(gcf, 'Figuren/opgave2j_mean_cv.png');

% Figuur 2: Variantie LOOCV-fout vs. N
figure('Position', [100 100 900 600]);
plot(N_values, var_cv, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Dataset-grootte N', 'FontSize', 12);
ylabel('Variantie van CV_{LOO}', 'FontSize', 12);
title(sprintf('Variantie LOOCV-fout vs. dataset-grootte (n=%d)', n), 'FontSize', 13);
grid on;
xlim([min(N_values)-10, max(N_values)+10]);

saveas(gcf, 'Figuren/opgave2j_var_cv.eps', 'epsc');
saveas(gcf, 'Figuren/opgave2j_var_cv.png');

% Samenvatting
fprintf('=== Samenvatting ===\n');
fprintf('N\t\tGem CV_LOO\tVar CV_LOO\n');
fprintf('---\t\t----------\t---------\n');
for i_N = 1:num_N
    fprintf('%d\t\t%.6f\t%.8f\n', N_values(i_N), mean_cv(i_N), var_cv(i_N));
end

% ========== HULPFUNCTIE ==========
function cv_loo = compute_loocv_subset(X, b, n)
    % Bereken LOOCV-fout voor gegeven dataset X, labels b, en modelgraad n
    N = size(X, 1);
    n_wrong = 0;
    
    % Precompute alle features
    Phi_all = build_poly_features(X(:,1), X(:,2), n);
    
    % Warmstart initialisatie
    beta_init = zeros(2*n + 1, 1);
    
    for i = 1:N
        % Laat punt i achterwege
        mask = true(N, 1);
        mask(i) = false;
        
        Phi_train = Phi_all(mask, :);
        b_train = b(mask);
        Phi_test = Phi_all(i, :);
        b_test = b(i);
        
        % Fit logistische regressie
        f = @(x) logreg_cost_pm1(x, Phi_train, b_train);
        df = @(x) logreg_grad_pm1(x, Phi_train, b_train);
        [beta, ~] = GD(f, df, beta_init, 1, 1e-5, 4000, 1e-4, 0.5);
        beta_init = beta; % Warmstart
        
        % Voorspel op testpunt
        b_pred = sign(Phi_test * beta);
        if b_pred == 0, b_pred = 1; end
        
        if b_pred ~= b_test
            n_wrong = n_wrong + 1;
        end
    end
    
    cv_loo = n_wrong / N;
end
