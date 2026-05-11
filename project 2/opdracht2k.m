% OPDRACHT 2K
% Bias van K-fold CV vs. aantal vouwen K
% Vaste parameters: N=400, n=3; variabel: K=1..15
% Voor elke K: 200 shuffles van de volledige dataset

clear; close all;

% Vaste parameters
N = 400;
n = 3;   % Modelgraad - gemakkelijk aan te passen
K_values = 1:15;
num_K = length(K_values);
num_samples = 200;

% Laad de dataset
load('DatasetCV.mat', 'x', 'y', 'cat');
x1 = x(:);
x2 = y(:);
b = cat(:);
X = [x1, x2];

if size(X, 1) ~= N
    error('DatasetCV.mat moet exact %d samples bevatten', N);
end

mean_cv = zeros(1, num_K);
var_cv = zeros(1, num_K);

tic;
fprintf('Opdracht 2(k): %d K-waarden, %d samples per K, N=%d, n=%d\n', num_K, num_samples, N, n);

for i_K = 1:num_K
    K = K_values(i_K);
    fprintf('K = %2d (%2d/%2d): ', K, i_K, num_K);
    
    cv_values = zeros(1, num_samples);
    
    % Genereer 200 shuffles van de volledige dataset
    for i_sample = 1:num_samples
        idx_shuffle = randperm(N);
        X_shuffled = X(idx_shuffle, :);
        b_shuffled = b(idx_shuffle);
        
        % Bereken K-fold CV op deze shuffle
        cv_k = compute_kfold_cv(X_shuffled, b_shuffled, K, n);
        cv_values(i_sample) = cv_k;
        
        if mod(i_sample, 50) == 0
            fprintf('.');
        end
    end
    
    mean_cv(i_K) = mean(cv_values);
    var_cv(i_K) = var(cv_values);
    
    fprintf(' mean=%.4f, var=%.6f\n', mean_cv(i_K), var_cv(i_K));
end

t = toc;
fprintf('Totale tijd: %.1f seconden\n\n', t);

% Figuur: Gemiddelde CV_K vs. K
figure('Position', [100 100 900 600]);
plot(K_values, mean_cv, 'b-o', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('Aantal vouwen K', 'FontSize', 12);
ylabel('Gemiddelde CV_K', 'FontSize', 12);
title(sprintf('Gemiddelde K-fold CV-fout vs. aantal vouwen (N=%d, n=%d)', N, n), 'FontSize', 13);
grid on;
xlim([0.5, 15.5]);
set(gca, 'XTick', K_values);

saveas(gcf, 'Figuren/opgave2k_mean_cv.eps', 'epsc');
saveas(gcf, 'Figuren/opgave2k_mean_cv.png');

% Samenvatting
fprintf('=== Samenvatting ===\n');
fprintf('K\tGem CV_K\tVar CV_K\n');
fprintf('---\t--------\t--------\n');
for i_K = 1:num_K
    fprintf('%d\t%.6f\t%.8f\n', K_values(i_K), mean_cv(i_K), var_cv(i_K));
end

% ========== HULPFUNCTIE ==========
function cv_k = compute_kfold_cv(X, b, K, n)
    % Bereken K-fold CV-fout voor gegeven dataset, labels, aantal vouwen en modelgraad
    N = size(X, 1);
    
    % Speciaal geval: K=1 betekent enkele random 50/50 split
    if K == 1
        split_idx = floor(N/2);
        X_train = X(1:split_idx, :);
        b_train = b(1:split_idx);
        X_test = X(split_idx+1:end, :);
        b_test = b(split_idx+1:end);
        
        Phi_train = build_poly_features(X_train(:,1), X_train(:,2), n);
        Phi_test = build_poly_features(X_test(:,1), X_test(:,2), n);
        
        f = @(x) logreg_cost_pm1(x, Phi_train, b_train);
        df = @(x) logreg_grad_pm1(x, Phi_train, b_train);
        [beta, ~] = GD(f, df, zeros(2*n+1, 1), 1, 1e-6, 15000, 1e-4, 0.5);
        
        n_wrong = 0;
        for i = 1:length(b_test)
            b_pred = sign(Phi_test(i, :) * beta);
            if b_pred == 0, b_pred = 1; end
            if b_pred ~= b_test(i)
                n_wrong = n_wrong + 1;
            end
        end
        
        cv_k = n_wrong / length(b_test);
        
    else
        % Standaard K-fold CV
        fold_size = floor(N / K);
        n_wrong = 0;
        total_test = 0;
        
        for k = 1:K
            % Bepaal grens van huidige vouw
            if k < K
                test_idx = (k-1)*fold_size + 1 : k*fold_size;
            else
                test_idx = (k-1)*fold_size + 1 : N;
            end
            
            % Train/test split
            mask = true(N, 1);
            mask(test_idx) = false;
            
            X_train = X(mask, :);
            b_train = b(mask);
            X_test = X(test_idx, :);
            b_test = b(test_idx);
            
            % Fit model
            Phi_train = build_poly_features(X_train(:,1), X_train(:,2), n);
            Phi_test = build_poly_features(X_test(:,1), X_test(:,2), n);
            
            f = @(x) logreg_cost_pm1(x, Phi_train, b_train);
            df = @(x) logreg_grad_pm1(x, Phi_train, b_train);
            [beta, ~] = GD(f, df, zeros(2*n+1, 1), 1, 1e-6, 15000, 1e-4, 0.5);
            
            % Tel missings
            for i = 1:length(b_test)
                b_pred = sign(Phi_test(i, :) * beta);
                if b_pred == 0, b_pred = 1; end
                if b_pred ~= b_test(i)
                    n_wrong = n_wrong + 1;
                end
            end
            
            total_test = total_test + length(b_test);
        end
        
        cv_k = n_wrong / total_test;
    end
end
