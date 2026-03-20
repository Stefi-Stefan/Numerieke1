% Opgave 5
% Benadering van twee datasets op hetzelfde rooster met n = m = 7.
clear; clc; close all;

% Gevraagde instellingen.
N = 31;
M = 31;
n = 7;
m = 7;

x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X,Y] = meshgrid(x,y);

%% Dataset 1: exp(1 - (2.5x)^2 - (2y)^2)
F1 = exp(1 - (2.5*X).^2 - (2*Y).^2);

C1 = kkb(x, y, F1, n, m);
A = get_leg_mtx(x, n);
B = get_leg_mtx(y, m);
Z1 = B * C1 * A';

figure;
scatter3(X(:), Y(:), F1(:), 12, 'k', 'filled'); hold on;
surf(X, Y, Z1, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
colormap turbo;
grid on;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Dataset 1: datapunten + veeltermbenadering (n=m=7)');
legend('Datapunten', 'Benaderend oppervlak', 'Location', 'best');
view(45,30);

%% Dataset 2: membrane(1,15)
F2 = membrane(1,15);

C2 = kkb(x, y, F2, n, m);
Z2 = B * C2 * A';

figure;
scatter3(X(:), Y(:), F2(:), 12, 'k', 'filled'); hold on;
surf(X, Y, Z2, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
colormap turbo;
grid on;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Dataset 2 (membrane): datapunten + veeltermbenadering (n=m=7)');
legend('Datapunten', 'Benaderend oppervlak', 'Location', 'best');
view(45,30);

%% Korte output met foutmaten
err1 = norm(F1 - Z1, 'fro');
err2 = norm(F2 - Z2, 'fro');

fprintf('Opgave 5 afgerond met n=m=7\n');
fprintf('||F1 - Z1||_F = %.3e\n', err1);
fprintf('||F2 - Z2||_F = %.3e\n', err2);
