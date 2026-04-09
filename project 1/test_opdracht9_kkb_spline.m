% testscript voor oefening 9
clear; clc; close all;

rng(1); % vaste seed voor reproduceerbaarheid

% we kiezen een kubische splinebasis op [0,1]
k = 3;
t = [0;0;0;0;0.15;0.35;0.55;0.75;1;1;1;1];
nb = length(t) - k - 1;
xplot = linspace(0,1,500);

%% TEST 1a: exacte reconstructie zonder ruis
% maak data van een gekende spline en controleer of die teruggevonden wordt

c_true = [1.2; -0.8; 0.3; 1.7; -1.1; 0.6; 0.2; -0.4];

R1 = 80;
x1 = linspace(0,1,R1)';
f1 = deboor(k, t, c_true, x1);
f1 = f1(:);

Z1 = kkb_spline(t, x1, f1, xplot);
Ztrue = deboor(k, t, c_true, xplot);

errInf_test1 = max(abs(Z1 - Ztrue));
err2_test1 = norm(Z1 - Ztrue, 2);

fprintf('\nTEST 1 (exacte reconstructie)\n');
fprintf('  ||Z - Ztrue||_inf = %.3e\n', errInf_test1);
fprintf('  ||Z - Ztrue||_2   = %.3e\n', err2_test1);

figure;
plot(xplot, Ztrue, 'k-', 'LineWidth', 1.7); hold on;
plot(xplot, Z1, 'r--', 'LineWidth', 1.3);
plot(x1, f1, 'bo', 'MarkerSize', 3);
grid on;
title('Test 1a: exacte reconstructie zonder ruis');
legend('ware spline','kkb\_spline','meetpunten','Location','best');
xlabel('x'); ylabel('y');

%% TEST 1b: residu op de trainingspunten
% vergelijk fitwaarden in de meetpunten met de data f1

Z1_train = kkb_spline(t, x1, f1, x1');
res_train_inf = max(abs(Z1_train(:) - f1));
res_train_2 = norm(Z1_train(:) - f1, 2);
res_train_rmse = sqrt(mean((Z1_train(:) - f1).^2));

fprintf('\nTEST 1b (residu op trainingspunten)\n');
fprintf('  ||residu||_inf = %.3e\n', res_train_inf);
fprintf('  ||residu||_2   = %.3e\n', res_train_2);
fprintf('  RMSE residu    = %.3e\n', res_train_rmse);

%% TEST 2: ruis op de meetdata
% voeg ruis toe en vergelijk de fit met de onderliggende ware spline

sigma = 0.05;
f2 = f1 + sigma*randn(size(f1));
Z2 = kkb_spline(t, x1, f2, xplot);

rmse_noisy_data = sqrt(mean((f2 - f1).^2));
rmse_fit_to_true = sqrt(mean((Z2 - Ztrue).^2));

fprintf('\nTEST 2 (ruisrobuustheid)\n');
fprintf('  RMSE van ruwe meetdata t.o.v. waarheid = %.3e\n', rmse_noisy_data);
fprintf('  RMSE van splinefit t.o.v. waarheid     = %.3e\n', rmse_fit_to_true);

figure;
subplot(1,2,1);
plot(xplot, Ztrue, 'k-', 'LineWidth', 1.7); hold on;
plot(xplot, Z2, 'm-', 'LineWidth', 1.4);
plot(x1, f2, '.', 'Color', [0.7 0.7 0.7]);
axis square;
grid on;
title('Fit bij ruis op de meetdata');
legend('ware spline','kkb\_spline met ruisdata','ruisdata','Location','northeast');
xlabel('x'); ylabel('y');

subplot(1,2,2);
plot(xplot, abs(Z2 - Ztrue), 'r-', 'LineWidth', 1.5);
axis square;
grid on;
title('Absolute benaderingsfout bij ruis');
legend('|Z_{fit} - Z_{waar}|', 'Location','best');
xlabel('x'); ylabel('Absolute fout');

%% TEST extra: invloed van het aantal meetpunten R
% fout daalt als er meer meetpunten zijn, en is perfect op machineprecisie
% na vanaf 8

Rvals = [3 4 5 6 7 8 9 10 20 40 80 160 320];
errInf_R = zeros(size(Rvals));
err2_R = zeros(size(Rvals));

for p = 1:length(Rvals)
    Rp = Rvals(p);
    xp = linspace(0,1,Rp)';
    fp = deboor(k, t, c_true, xp);
    fp = fp(:);

    Zp = kkb_spline(t, xp, fp, xplot);
    errInf_R(p) = max(abs(Zp - Ztrue));
    err2_R(p) = norm(Zp - Ztrue, 2);
end

fprintf('\nTEST extra (invloed aantal meetpunten)\n');
fprintf('  R        ||Z-Ztrue||_inf        ||Z-Ztrue||_2\n');
for p = 1:length(Rvals)
    fprintf('  %-4d     %-16.3e     %-16.3e\n', Rvals(p), errInf_R(p), err2_R(p));
end

figure;
loglog(Rvals, errInf_R, 'o-', 'LineWidth', 1.5); hold on;
loglog(Rvals, err2_R, 's-', 'LineWidth', 1.5);
grid on;
title('Test extra: fout versus aantal meetpunten');
legend('||Z-Ztrue||_\infty','||Z-Ztrue||_2','Location','eastoutside');
xlabel('R');
ylabel('fout');

%% TEST 3: randgevallen - welke functies werken goed/slecht?
% we gebruiken telkens dezelfde knopenrij t en hetzelfde aantal meetpunten
% zo vergelijken we eerlijk hoe "moeilijk" de functie is voor de splinebasis

R3 = 120;
x3 = linspace(0,1,R3)';

% 1) zeer glad (goed benaderbaar)
f_smooth = exp(x3) .* cos(2*pi*x3);

% 2) knik in x=0.5 (minder glad)
f_kink = abs(x3 - 0.5);

% 3) hoge frequentie (moeilijk met beperkt aantal knopen)
f_osc = sin(20*pi*x3);

% 4) sprongfunctie (zeer lastig voor gladde splines)
f_step = double(x3 >= 0.5);

% evaluatie op fijn rooster voor foutmeting en figuren
f_smooth_plot = exp(xplot) .* cos(2*pi*xplot);
f_kink_plot = abs(xplot - 0.5);
f_osc_plot = sin(20*pi*xplot);
f_step_plot = double(xplot >= 0.5);

Z_smooth = kkb_spline(t, x3, f_smooth, xplot);
Z_kink = kkb_spline(t, x3, f_kink, xplot);
Z_osc = kkb_spline(t, x3, f_osc, xplot);
Z_step = kkb_spline(t, x3, f_step, xplot);

rmse_smooth = sqrt(mean((Z_smooth - f_smooth_plot).^2));
rmse_kink = sqrt(mean((Z_kink - f_kink_plot).^2));
rmse_osc = sqrt(mean((Z_osc - f_osc_plot).^2));
rmse_step = sqrt(mean((Z_step - f_step_plot).^2));

max_smooth = max(abs(Z_smooth - f_smooth_plot));
max_kink = max(abs(Z_kink - f_kink_plot));
max_osc = max(abs(Z_osc - f_osc_plot));
max_step = max(abs(Z_step - f_step_plot));

fprintf('\nTEST 3 (randgevallen: functie-afhankelijk gedrag)\n');
fprintf('  Functie                     RMSE                maxfout\n');
fprintf('  glad exp(x)cos(2pix)       %-18.3e  %-18.3e\n', rmse_smooth, max_smooth);
fprintf('  knik |x-0.5|               %-18.3e  %-18.3e\n', rmse_kink, max_kink);
fprintf('  hoogfreq sin(20pix)        %-18.3e  %-18.3e\n', rmse_osc, max_osc);
fprintf('  sprong 1_{x>=0.5}          %-18.3e  %-18.3e\n', rmse_step, max_step);

figure;
subplot(2,2,1);
plot(xplot, f_smooth_plot, 'k-', 'LineWidth', 1.5); hold on;
plot(xplot, Z_smooth, 'r--', 'LineWidth', 1.3); grid on;
title('Glad: exp(x)cos(2\pix)');
legend('waar','fit','Location','north');
xlabel('x'); ylabel('y');

subplot(2,2,2);
plot(xplot, f_kink_plot, 'k-', 'LineWidth', 1.5); hold on;
plot(xplot, Z_kink, 'r--', 'LineWidth', 1.3); grid on;
title('Knik: |x-0.5|');
legend('waar','fit','Location','north');
xlabel('x'); ylabel('y');

subplot(2,2,3);
plot(xplot, f_osc_plot, 'k-', 'LineWidth', 1.5); hold on;
plot(xplot, Z_osc, 'r--', 'LineWidth', 1.3); grid on;
title('Hoge frequentie: sin(20\pix)');
legend('waar','fit','Location','north');
xlabel('x'); ylabel('y');

subplot(2,2,4);
plot(xplot, f_step_plot, 'k-', 'LineWidth', 1.5); hold on;
plot(xplot, Z_step, 'r--', 'LineWidth', 1.3); grid on;
title('Sprongfunctie');
legend('waar','fit','Location','north');
xlabel('x'); ylabel('y');
%% TEST 4: randpunten x=0 en x=1
% deze punten geven vaak problemen bij splines met herhaalde eindknopen

x_edge = [0 1];
f_edge_true = deboor(k, t, c_true, x_edge);
f_edge_fit = kkb_spline(t, x1, f1, x_edge);
err_edge = abs(f_edge_fit - f_edge_true);

fprintf('\nTEST 4 (randpunten)\n');
fprintf('  x = 0: |f_fit - f_true| = %.3e\n', err_edge(1));
fprintf('  x = 1: |f_fit - f_true| = %.3e\n', err_edge(2));

%% Samenvattingstabel
fprintf('\nSAMENVATTING TESTEN\n');
fprintf('  Test                                Maat                        Waarde\n');
fprintf('  Test 1 exacte reconstructie         ||Z-Ztrue||_inf             %.3e\n', errInf_test1);
fprintf('  Test 1 exacte reconstructie         ||Z-Ztrue||_2               %.3e\n', err2_test1);
fprintf('  Test 1b trainingsresidu             RMSE                        %.3e\n', res_train_rmse);
fprintf('  Test 2 ruis                         RMSE fit tov waarheid       %.3e\n', rmse_fit_to_true);
fprintf('  Test 4 randpunten                   max randpuntfout            %.3e\n', max(err_edge));