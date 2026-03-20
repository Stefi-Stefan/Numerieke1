% Test voor opgave 4
% Correctheidstest met een voorbeeld waarvoor de oplossing gekend is.
clear; clc; close all;

% Kies rooster op [-1,1].
N = 11;
M = 9;
x = linspace(-1,1,N);
y = linspace(-1,1,M);

% Kies graden.
n = 3;
m = 2;

% Maak een gekende coefficientenmatrix C_true.
C_true = [ 1.0  -0.4   0.2   0.1;
          -0.3   0.7  -0.5   0.0;
           0.4   0.2   0.0  -0.2];

% Bouw F exact op uit de Legendrebasis: F = B * C_true * A^T
A = get_leg_mtx(x, n);
B = get_leg_mtx(y, m);
F = B * C_true * A';

% Herstel C met de implementatie.
C_est = kkb(x, y, F, n, m);

% Vergelijk met exacte C.
errInf = max(abs(C_est(:) - C_true(:)));
errFro = norm(C_est - C_true, 'fro');

fprintf('Test opgave 4 (gekende oplossing)\n');
fprintf('||C_est - C_true||_inf = %.3e\n', errInf);
fprintf('||C_est - C_true||_F   = %.3e\n', errFro);

% Extra check op reconstructie van F.
F_rec = B * C_est * A';
recErr = norm(F_rec - F, 'fro');
fprintf('||F_rec - F||_F        = %.3e\n', recErr);
