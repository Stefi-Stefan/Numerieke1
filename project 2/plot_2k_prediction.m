% plot_2k_prediction.m
% Plot mean and variance for K=1..15 using provided results
% with a predicted value for K=15.

clear; close all; clc;

K = 1:15;
mean_vals = [0.1174, 0.1169, 0.1158, 0.1160, 0.1164, 0.1159, 0.1168, 0.1165, 0.1167, 0.1167, 0.1169, 0.1164, 0.1172, 0.1169, 0.1168]; % predicted K=15 -> 0.1168
var_vals  = [0.000348, 0.000057, 0.000053, 0.000040, 0.000031, 0.000029, 0.000032, 0.000029, 0.000026, 0.000026, 0.000024, 0.000021, 0.000021, 0.000019, 0.000018]; % predicted K=15 -> 1.8e-05

std_vals = sqrt(var_vals);

% Create figure with two subplots
fig = figure('Color','w','Position',[100 100 1000 420]);

subplot(1,2,1);
errorbar(K, mean_vals, std_vals, 'o-','LineWidth',1.4,'MarkerSize',6);
xlabel('Aantal vouwen K');
ylabel('Gemiddelde CV-fout (mean)');
title('Gemiddelde CV-fout');
grid on;
xticks(K);

subplot(1,2,2);
plot(K, var_vals, 's-','LineWidth',1.4,'MarkerSize',6);
xlabel('Aantal vouwen K');
ylabel('Varianties van CV-fout');
title('Varianties per K (sample variance over shuffles)');
grid on;
xticks(K);

outDir = fullfile(pwd, 'Figuren');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

saveas(fig, fullfile(outDir, 'opgave2k_predicted.png'));
saveas(fig, fullfile(outDir, 'opgave2k_predicted.eps'), 'epsc');

fprintf('Voorspelde voor K=15: mean=%.4f, var=%.6f\n', mean_vals(end), var_vals(end));
