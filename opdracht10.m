% Opgave 10
% Vergelijking van twee knopenrijen t1 en t2 voor splinefit van meetgegevens.
clear; clc; close all;

% Laad data uit het meegegeven bestand.
S = load('exercise10.mat');

% Verwachte variabelen volgens opgave.
x = S.x;
f = S.y;
t1 = S.t1;
t2 = S.t2;

% Zorg voor correcte vorm.
x = x(:);
f = f(:);
t1 = t1(:);
t2 = t2(:);

% Fijn rooster voor visualisatie.
xplot = linspace(min(x), max(x), 600);

% Bereken beide splinebenaderingen.
Z1 = kkb_spline(t1, x, f, xplot);
Z2 = kkb_spline(t2, x, f, xplot);

% Evalueer ook op de meetpunten voor kwantitatieve vergelijking.
F1_train = kkb_spline(t1, x, f, x');
F2_train = kkb_spline(t2, x, f, x');

% Foutmaten op de meetpunten.
rmse1 = sqrt(mean((F1_train(:) - f).^2));
rmse2 = sqrt(mean((F2_train(:) - f).^2));
maxerr1 = max(abs(F1_train(:) - f));
maxerr2 = max(abs(F2_train(:) - f));

% Plot meetpunten met beide fits samen.
figure;
plot(x, f, 'ko', 'MarkerSize', 4, 'DisplayName', 'meetpunten'); hold on;
plot(xplot, Z1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'fit met t1');
plot(xplot, Z2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'fit met t2');

% Toon de knopen als kleine markeringen onderaan de figuur.
yl = ylim;
dy = yl(2) - yl(1);
t1_plot = t1(t1 >= min(xplot) & t1 <= max(xplot));
t2_plot = t2(t2 >= min(xplot) & t2 <= max(xplot));
plot(t1_plot, (yl(1) - 0.05*dy) * ones(size(t1_plot)), 'bv', 'MarkerSize', 4, 'HandleVisibility', 'off');
plot(t2_plot, (yl(1) - 0.10*dy) * ones(size(t2_plot)), 'rv', 'MarkerSize', 4, 'HandleVisibility', 'off');
ylim([yl(1) - 0.15*dy, yl(2)]);

grid on;
xlabel('x'); ylabel('f(x)');
title('Opgave 10: vergelijking van t1 en t2');
legend('Location', 'best');

% Extra figuur met aparte panelen.
figure;
subplot(2,1,1);
plot(x, f, 'ko', 'MarkerSize', 4); hold on;
plot(xplot, Z1, 'b-', 'LineWidth', 1.5);

yl = ylim;
dy = yl(2) - yl(1);
t1_plot = t1(t1 >= min(xplot) & t1 <= max(xplot));
plot(t1_plot, (yl(1) - 0.05*dy) * ones(size(t1_plot)), 'bv', 'MarkerSize', 4);
ylim([yl(1) - 0.10*dy, yl(2)]);

grid on;
title('Splinefit met knopenrij t1');
xlabel('x'); ylabel('f(x)');

subplot(2,1,2);
plot(x, f, 'ko', 'MarkerSize', 4); hold on;
plot(xplot, Z2, 'r--', 'LineWidth', 1.5);

yl = ylim;
dy = yl(2) - yl(1);
t2_plot = t2(t2 >= min(xplot) & t2 <= max(xplot));
plot(t2_plot, (yl(1) - 0.05*dy) * ones(size(t2_plot)), 'rv', 'MarkerSize', 4);
ylim([yl(1) - 0.10*dy, yl(2)]);

grid on;
title('Splinefit met knopenrij t2');
xlabel('x'); ylabel('f(x)');

% Toon kwantitatieve vergelijking in command window.
fprintf('\nKWANTITATIEVE VERGELIJKING OP MEETPUNTEN\n');
fprintf('t1: RMSE = %.3e, maxfout = %.3e\n', rmse1, maxerr1);
fprintf('t2: RMSE = %.3e, maxfout = %.3e\n', rmse2, maxerr2);

if rmse1 < rmse2
    fprintf('Beste benadering volgens RMSE: t1\n');
elseif rmse2 < rmse1
    fprintf('Beste benadering volgens RMSE: t2\n');
else
    fprintf('Volgens RMSE zijn t1 en t2 gelijkwaardig.\n');
end
