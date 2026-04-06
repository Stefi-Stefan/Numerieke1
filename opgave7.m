clear; clc; close all;

F = double(imread('matterhorn.png'));
[M, N] = size(F);
x = linspace(-1,1,N);
y = linspace(-1,1,M);
[X, Y] = meshgrid(x, y);
n = 25;
m = 25;
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

%% bovenaanzicht
figure;
subplot(1,2,1);
imagesc(x, y, F);
set(gca, 'YDir', 'normal');
axis square;
colorbar;
title('imagesc: originele data');
xlabel('x'); ylabel('y');

subplot(1,2,2);
imagesc(x, y, Z);
set(gca, 'YDir', 'normal');
axis square;
colorbar;
title('imagesc: veeltermbenadering');
xlabel('x'); ylabel('y');

%% contour
figure;
subplot(1,2,1);
contourf(X, Y, F,10, 'LineColor', 'none');
axis square;
colorbar;
title('contourf: originele data');
xlabel('x'); ylabel('y');

subplot(1,2,2);
contourf(X, Y, Z,10, 'LineColor', 'none');
axis square;
colorbar;
title('contourf: veeltermbenadering');
xlabel('x'); ylabel('y');

%% grootte fout
errFro = norm(F - Z, 'fro');
relErr = errFro / max(norm(F, 'fro'), eps);
fprintf('Opgave 7 benaderd met n=m=25\n');
fprintf('||F - Z||_F = %.3e\n', errFro);
fprintf('Relatieve fout = %.3e\n', relErr);

errAbs = abs(F - Z);
figure;
imagesc(x, y, flipud(errAbs));
set(gca, 'YDir', 'normal'); 
axis square;
colorbar;
colormap(hot); 
title('Absolute Fout abs(F - Z)');
xlabel('x'); ylabel('y');
