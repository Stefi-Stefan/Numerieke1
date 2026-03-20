% Opgave 7
% Veeltermbenadering van de Matterhorn-data met n = m = 25.
clear; clc; close all;

% Lees data in.
F_raw = imread('matterhorn.png');

% Zet om naar 2D matrix indien nodig (grayscale).
if ndims(F_raw) == 3
    F = double(rgb2gray(F_raw));
else
    F = double(F_raw);
end

% Afmetingen en genormaliseerd rooster.
[M, N] = size(F);
x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X, Y] = meshgrid(x, y);

% Gevraagde graden.
n = 25;
m = 25;

% Bepaal de benadering.
C = kkb(x, y, F, n, m);
A = get_leg_mtx(x, n);
B = get_leg_mtx(y, m);
Z = B * C * A';

%% 3D-oppervlakken: origineel en benadering
figure;
subplot(1,2,1);
surf(X, Y, F, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
xlim([-1,1]); ylim([-1,1]); zlim([0,250]);
title('Matterhorn - originele data');
xlabel('x'); ylabel('y'); zlabel('hoogte');
view(45,30);
camlight headlight;

subplot(1,2,2);
surf(X, Y, Z, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
xlim([-1,1]); ylim([-1,1]); zlim([0,250]);
title('Matterhorn - veeltermbenadering (n=m=25)');
xlabel('x'); ylabel('y'); zlabel('hoogte');
view(45,30);
camlight headlight;

%% Bovenaanzicht met imagesc
figure;
subplot(1,2,1);
imagesc(F);
axis image;
colorbar;
title('imagesc: originele data');

subplot(1,2,2);
imagesc(Z);
axis image;
colorbar;
title('imagesc: veeltermbenadering');

%% Contourplots
figure;
subplot(1,2,1);
contourf(X, Y, F, 25, 'LineColor', 'none');
colorbar;
title('contourf: originele data');
xlabel('x'); ylabel('y');

subplot(1,2,2);
contourf(X, Y, Z, 25, 'LineColor', 'none');
colorbar;
title('contourf: veeltermbenadering');
xlabel('x'); ylabel('y');

%% Foutmaat
errFro = norm(F - Z, 'fro');
relErr = errFro / max(norm(F, 'fro'), eps);

fprintf('Opgave 7 afgerond met n=m=25\n');
fprintf('||F - Z||_F = %.3e\n', errFro);
fprintf('Relatieve fout = %.3e\n', relErr);
