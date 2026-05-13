% OPDRACHT 2K - FAST VARIANT
% Versneld experiment voor bias van K-fold CV vs. K
% Wijzigingen t.o.v. origineel:
% - num_samples verlaagd (minder shuffles)
% - GD tolerantie verhoogd (minder precieze oplossing)
% - GD max_iter verlaagd (minder iteraties)
% Doel: veel lagere rekentijd, behoud van dezelfde trend.

clear; close all;

% Vaste parameters
N = 400;
n = 3;   % Modelgraad - gemakkelijk aan te passen
K_values = 1:15;
num_K = length(K_values);
num_samples = 80; % verlaagd van 200 naar 80 voor snellere runs

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

% Snellere GD-instellingen
gd_tol = 1e-5;      % ruimere tolerantie (minder iteraties)
gd_max_iter = 4000; % beperk aantal iteraties

fprintf('Opdracht 2(k) fast: %d K-waarden, %d samples per K, N=%d, n=%d\n', num_K, num_samples, N, n);

for i_K = 1:num_K
    K = K_values(i_K);
    fprintf('K = %2d (%2d/%2d): ', K, i_K, num_K);
    
    cv_values = zeros(1, num_samples);
    
    % Genereer shuffles van de volledige dataset
    for i_sample = 1:num_samples
        idx_shuffle = randperm(N);
        X_shuffled = X(idx_shuffle, :);
        b_shuffled = b(idx_shuffle);
        
        % Bereken K-fold CV op deze shuffle (met snellere GD instellingen)
        cv_k = compute_kfold_cv_fast(X_shuffled, b_shuffled, K, n, gd_tol, gd_max_iter);
        cv_values(i_sample) = cv_k;
        
        if mod(i_sample, 20) == 0
            fprintf('.');
        end
    end
    
    mean_cv(i_K) = mean(cv_values);
    var_cv(i_K) = var(cv_values);
    
    fprintf(' mean=%.4f, var=%.6f\n', mean_cv(i_K), var_cv(i_K));
end

% Figuur: Gemiddelde CV_K vs. K en variantie
figure('Position',[100 100 1000 420]);
subplot(1,2,1);
errorbar(K_values, mean_cv, sqrt(var_cv), 'b-o','LineWidth',1.6,'MarkerSize',6);
xlabel('Aantal vouwen K');
ylabel('Gemiddelde CV_K');
title(sprintf('Gemiddelde K-fold CV-fout vs. K (fast) -- n=%d', n));
grid on; xticks(K_values);

subplot(1,2,2);
plot(K_values, var_cv, 'r-s','LineWidth',1.6,'MarkerSize',6);
xlabel('Aantal vouwen K');
ylabel('Varianties van CV_K');
title('Varianties over shuffles');
grid on; xticks(K_values);

outDir = fullfile(pwd, 'Figuren');
if ~exist(outDir,'dir'), mkdir(outDir); end
saveas(gcf, fullfile(outDir, 'opgave2k_fast_mean_var.png'));

fprintf('\nKlaar: fast variant gemaakt. Results opgeslagen in Figuren/opgave2k_fast_mean_var.png\n');

% ===== HULPFUNCTIE (snellere GD-instellingen) =====
function cv_k = compute_kfold_cv_fast(X, b, K, n, gd_tol, gd_max_iter)
    N = size(X,1);
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
        [beta, ~] = GD(f, df, zeros(2*n+1,1), 1, gd_tol, gd_max_iter, 1e-4, 0.5);
        
        ypred = sign(Phi_test * beta);
        ypred(ypred==0) = 1;
        cv_k = sum(ypred ~= b_test) / numel(b_test);
    else
        fold_size = floor(N / K);
        n_wrong = 0;
        total_test = 0;
        for k = 1:K
            if k < K
                test_idx = (k-1)*fold_size + 1 : k*fold_size;
            else
                test_idx = (k-1)*fold_size + 1 : N;
            end
            mask = true(N,1);
            mask(test_idx) = false;
            X_train = X(mask,:); b_train = b(mask);
            X_test = X(test_idx,:); b_test = b(test_idx);
            
            Phi_train = build_poly_features(X_train(:,1), X_train(:,2), n);
            Phi_test = build_poly_features(X_test(:,1), X_test(:,2), n);
            
            f = @(x) logreg_cost_pm1(x, Phi_train, b_train);
            df = @(x) logreg_grad_pm1(x, Phi_train, b_train);
            [beta, ~] = GD(f, df, zeros(2*n+1,1), 1, gd_tol, gd_max_iter, 1e-4, 0.5);
            
            ypred = sign(Phi_test * beta);
            ypred(ypred==0) = 1;
            n_wrong = n_wrong + sum(ypred ~= b_test);
            total_test = total_test + numel(b_test);
        end
        cv_k = n_wrong / total_test;
    end
end
