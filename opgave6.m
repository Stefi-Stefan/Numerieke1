% Opgave 6
% Verloop van de kostfunctie voor n = m = 1,...,20 voor twee datasets.
clear; clc; close all;

% Zelfde rooster als opgave 5.
N = 31;
M = 31;
x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X,Y] = meshgrid(x,y);

% Datasets.
F1 = exp(1 - (2.5*X).^2 - (2*Y).^2);
F2 = membrane(1,15);

% Graden n = m = 1,...,20.
deg = 1:20;
J1 = zeros(size(deg));
J2 = zeros(size(deg));

for idx = 1:length(deg)
    d = deg(idx);

    C1 = kkb(x, y, F1, d, d);
    A = get_leg_mtx(x, d);
    B = get_leg_mtx(y, d);
    Z1 = B * C1 * A';

    C2 = kkb(x, y, F2, d, d);
    Z2 = B * C2 * A';

    % Kostfunctie J = sum_{i,j} (F - Z)^2
    E1 = F1 - Z1;
    E2 = F2 - Z2;
    J1(idx) = sum(E1(:).^2);
    J2(idx) = sum(E2(:).^2);
end

% Plot van de kostfuncties op dezelfde figuur.
figure;
semilogy(deg, J1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(deg, J2, 's-', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
xlabel('graad d (met n = m = d)');
ylabel('kostfunctie J(d) = \Sigma_{i,j}(F_{ij} - Z_{ij})^2');
title('Opgave 6: verloop van de kostfunctie');
legend('Dataset 1: exp(1-(2.5x)^2-(2y)^2)', 'Dataset 2: membrane(1,15)', 'Location', 'best');

% Kwantitatieve output voor bespreking.
ratio1 = J1(end) / J1(1);
ratio2 = J2(end) / J2(1);

fprintf('Opgave 6: kostfunctie voor n = m = 1,...,20\n');
fprintf('Dataset 1: J(1)=%.3e, J(20)=%.3e, verhouding J(20)/J(1)=%.3e\n', J1(1), J1(end), ratio1);
fprintf('Dataset 2: J(1)=%.3e, J(20)=%.3e, verhouding J(20)/J(1)=%.3e\n', J2(1), J2(end), ratio2);

% Extra indicatie van afneming in laatste stappen.
fprintf('Dataset 1: relatieve daling van d=19 naar d=20: %.3e\n', (J1(end-1)-J1(end))/max(J1(end-1),eps));
fprintf('Dataset 2: relatieve daling van d=19 naar d=20: %.3e\n', (J2(end-1)-J2(end))/max(J2(end-1),eps));
