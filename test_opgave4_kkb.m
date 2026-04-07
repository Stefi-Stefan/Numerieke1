%Correctheidstest voor opgave 4 met een voorbeeld waarvoor de oplossing gekend is
clear; clc; close all;

% bv 11x9 rooster op [-1,1].
N = 11;
M = 9;
x = linspace(-1,1,N);
y = linspace(-1,1,M);

n = 3; % gr x
m = 2; % gr y

% coefficientenmatrix C_true  (coefficienten vd veelterm in de
% Legendrebasis)
C_true = [ 1.0  -0.4   0.2   0.1;
    -0.3   0.7  -0.5   0.0;
    0.4   0.2   0.0  -0.2];

% F uit de Legendrebasis met diezelfde coeff : F = B * C_true * A^T
A = get_leg_mtx(x, n);
B = get_leg_mtx(y, m);
F = B * C_true * A';

C_est = kkb(x, y, F, n, m); % de C met kkb berekenen voor onze situatie


% deze berekende C_est zou bij benadering hetzelfde moeten zijn als de
% oorspronkelijk gekozen C_true
errInf = max(abs(C_est(:) - C_true(:)));
errFro = norm(C_est - C_true, 'fro');

fprintf('Test opgave 4 (gekende oplossing)\n');
fprintf('||C_est - C_true||_inf = %.3e\n', errInf);
fprintf('||C_est - C_true||_F   = %.3e\n', errFro);

% reconstructie van F zou ook bij benadering overeen moeten komen met de
% oorspronkelijk gemaakte data
F_rec = B * C_est * A';
recErr = norm(F_rec - F, 'fro');
fprintf('||F_rec - F||_F        = %.3e\n', recErr);
