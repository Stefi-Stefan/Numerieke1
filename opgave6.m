clear; clc; close all;

N = 31;
M = 31;
x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X,Y] = meshgrid(x,y);

F1 = exp(1 - (2.5*X).^2 - (2*Y).^2);
F2 = membrane(1,15);

degrees = 1:20; %graad n = m = 1,...,20.
K1 = zeros(size(degrees));
K2 = zeros(size(degrees));

for idx = 1:length(degrees)
    degree = degrees(idx);

    C1 = kkb(x, y, F1, degree, degree);
    A = get_leg_mtx(x, degree);
    B = get_leg_mtx(y, degree);
    Z1 = B * C1 * A';

    C2 = kkb(x, y, F2, degree, degree);
    Z2 = B * C2 * A';

    % kostfunctie K = sum_{i,j} (F_{ji} - Z_{ji})^2
    % = |F-Z|_F^2 = |vec(F)- vec(Z)|_2^2 = sum_{i,j} (f_{ij} - z(x_i,y_j})^2
    E1 = F1 - Z1;
    E2 = F2 - Z2;
    K1(idx) = sum(E1(:).^2);
    K2(idx) = sum(E2(:).^2);
end

figure;
semilogy(degrees, K1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(degrees, K2, 's-', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
xlabel('graad d (met n = m = d)');
ylabel('kostfunctie K(d) = \Sigma_{i,j} (f_{ij} - z(x_i,y_j))^2');
title('Opgave 6: verloop van de kostfunctie');
legend('Dataset 1: exp(1-(2.5x)^2-(2y)^2)', 'Dataset 2: membrane(1,15)', 'Location', 'northeast');

ratio1 = K1(end) / K1(1);
ratio2 = K2(end) / K2(1);

fprintf('Opgave 6: kostfunctie voor n = m = 1,...,20\n');
fprintf('Dataset 1: K(1)=%.3e, K(20)=%.3e, verhouding K(20)/K(1)=%.3e\n', K1(1), K1(end), ratio1);
fprintf('Dataset 2: K(1)=%.3e, K(20)=%.3e, verhouding K(20)/K(1)=%.3e\n', K2(1), K2(end), ratio2);
fprintf('Dataset 1: relatieve daling van d=19 naar d=20: %.3e\n', (K1(end-1)-K1(end))/max(K1(end-1),eps));
fprintf('Dataset 2: relatieve daling van d=19 naar d=20: %.3e\n', (K2(end-1)-K2(end))/max(K2(end-1),eps));
    