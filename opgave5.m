clear; clc; close all;

N = 31;
M = 31;
n = 7;
m = 7;

x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X,Y] = meshgrid(x,y);

%% Dataset 1: f_ij = exp(1 - (2.5*x_i)^2 - (2*y_j)^2)
F1 = exp(1 - (2.5*X).^2 - (2*Y).^2);

C1 = kkb(x, y, F1, n, m);
A = get_leg_mtx(x, n);
B = get_leg_mtx(y, m);
Z1 = B * C1 * A'; %in opgave stond hint over kron, wij hebben geen kron
% gebruikt maar wel de manier van Z = B * C * A'  die we in opgave 1
% bewezen hebben in het verslag omdat dit is veel efficienter is. met de 
% kron aanpak zou het worden:
%Gamma = kron(A, B);
%c1_vec = C1(:);
%z1_vec = Gamma * c1_vec;
%Z1 = reshape(z1_vec, M, N); %vector naar matrix terug

figure;
scatter3(X(:), Y(:), F1(:), 12, 'k', 'filled'); hold on;
surf(X, Y, Z1, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
colormap turbo;
grid on;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Dataset 1: datapunten + veeltermbenadering (n=m=7)');
legend('Datapunten', 'Benaderend oppervlak', 'Location', 'best');
view(45,30);

%% Dataset 2: F = membrane(1,15)
F2 = membrane(1,15);

C2 = kkb(x, y, F2, n, m);
Z2 = B * C2 * A';

figure;
scatter3(X(:), Y(:), F2(:), 12, 'k', 'filled'); hold on;
surf(X, Y, Z2, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
colormap turbo;
grid on;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Dataset 2 (membrane): datapunten + veeltermbenadering (n=m=7)');
legend('Datapunten', 'Benaderend oppervlak', 'Location', 'best');
view(-160,20);

%% grootte fouten
err1 = norm(F1 - Z1, 'fro');
err2 = norm(F2 - Z2, 'fro');

fprintf('Opgave 5 afgerond met n=m=7\n');
fprintf('||F1 - Z1||_F = %.3e\n', err1);
fprintf('||F2 - Z2||_F = %.3e\n', err2);
