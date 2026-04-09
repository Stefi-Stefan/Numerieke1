clear; clc; close all;

%% DEELVRAAG 1 Gram-Schmidt: Algoritmes
% Zie herGS.m, klGS.m en modGS.m


%% DEELVRAAG 2 Gram-Schmidt: Equivalentie
% Zie verslag


%% DEELVRAAG 3 Gram-Schmidt: Matrix A genereren
m = 200;
n = 40;
Q = orth(randn(m,n)); % Random orthonormale matrix (m x n)
D = diag(2.^-(0:n-1)); % Diagonaalmatrix met D(i,i) = 2^(1-i)
V = eye(n);
V(1, :) = 1; % Matrix met op de diagonaal en de eerste rij enen 
A = Q * D * V;


%% DEELVRAAG 4 Gram-Schmidt: Toepassen van Klassieke en Gewijzigde Gram-Schmidt
% !!! Opmerking: Om deze vraag te kunnen runnen, moet eerst sectie deelvraag3
% gerund worden !!!
j_waarde = 1:n;
min_y = 1e-16; % Minimum y-waarde voor de grafiek
max_y = 1e1; % Maximum y-waarde voor de grafiek

% Toepassen van de algoritmes op de matrix A
[Q_kl, ~] = klGS(A); % Klassieke GS
[Q_mod, ~] = modGS(A); % Gewijzigde GS

% Berekenen van de foutnorm per kolom
error_kl = vecnorm(Q_kl - Q, 2, 1); % ||q̃j - qj||2 voor klassieke GS
error_mod = vecnorm(Q_mod - Q, 2, 1); % ||q̃j - qj||2 voor gewijzigde GS

% Plotten van de foutnormen
figure;
semilogy(j_waarde, error_kl, '.-', 'LineWidth', 1.5, 'DisplayName', 'Klassieke GS - Fout');
hold on;
semilogy(j_waarde, error_mod, '.-', 'LineWidth', 1.5, 'DisplayName', 'Gewijzigde GS - Fout');
xlabel('Kolomindex $j$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Foutnorm $\|\tilde{q}_j - q_j\|_2$', 'Interpreter', 'latex', 'FontSize', 12);
ylim([min_y max_y]);
title('\textbf{Foutvergelijking tussen Klassieke en Gewijzigde GS}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 10, 'Interpreter', 'latex');
grid on;

%% Deelvraag 5 Gram-Schmidt: Foutenanalyse van Klassieke en Gewijzigde Gram-Schmidt

% Parameters
m = 200; % Aantal rijen van A
n = 40;  % Aantal kolommen van A
epsmach = eps; % Machineprecisie

% Genereren van de matrix A 
Q = orth(randn(m, n)); % Random orthonormale matrix (m x n)
D = diag(2.^-(0:n-1)); % Diagonaalmatrix met D(i,i) = 2^(1-i)
V = eye(n);
V(1, :) = 1; % Eerste rij vol met enen
A = Q * D * V;

Q_orig = Q; % Bewaar originele Q voor vergelijking

% Toepassen van Klassieke en Gewijzigde Gram-Schmidt
[Q_kl, ~] = klGS(A);    % Klassieke GS
[Q_mod, ~] = modGS(A);  % Gewijzigde GS

% Berekenen van de werkelijke fouten (verschil tussen Q en Q_orig)
error_kl = vecnorm(Q_kl - Q_orig, 2, 1);   % ||q̃_j - q_j||_2 voor klassieke GS
error_mod = vecnorm(Q_mod - Q_orig, 2, 1); % ||q̃_j - q_j||_2 voor gewijzigde GS

% Definieer theoretische foute als een anonieme functies
theor_fout_mod_func = @(i, Q_th_mod, epsmach) vecnorm(-2.^(i-1) * epsmach .* Q_th_mod(:,1), 2, 1);
theor_fout_kl_func = @(i, Q_th_kl, epsmach) vecnorm(-2.^(i-1) * epsmach .* Q_th_kl(:,1) + sum(Q_th_kl(:,2:i-1) .* (2.^(i + (2:i-1) - 2) * epsmach), 2), 2, 1);

% Theoretische fouten berekenen
theor_fout_kl = arrayfun(@(i) theor_fout_kl_func(i, Q_orig, epsmach), 1:n);
theor_fout_mod = arrayfun(@(i) theor_fout_mod_func(i, Q_orig, epsmach), 1:n);

% Plotten van de werkelijke fout en theoretische fout voor Klassieke GS
figure;
semilogy(1:n, error_kl, '.-', 'LineWidth', 1.5, 'DisplayName', 'Experimentele Fout');
hold on;
semilogy(1:n, theor_fout_kl, '--k', 'LineWidth', 1.5, 'DisplayName', 'Theoretische fout');
xlabel('Kolomindex $j$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Foutnorm $\|\tilde{q}_j - q_j\|_2$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Foutenanalyse: Klassieke Gram-Schmidt}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 10, 'Interpreter', 'latex');
grid on;

% Plotten van de werkelijke fout en theoretische fout voor Gewijzigde GS
figure;
semilogy(1:n, error_mod, '.-r', 'LineWidth', 1.5, 'DisplayName', 'Experimentele Fout');
hold on;
semilogy(1:n, theor_fout_mod, '--k', 'LineWidth', 1.5, 'DisplayName', 'Theoretische fout');
xlabel('Kolomindex $j$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Foutnorm $\|\tilde{q}_j - q_j\|_2$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Foutenanalyse: Gewijzigde Gram-Schmidt}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 10, 'Interpreter', 'latex');
grid on;


%% DEELVRAAG 6 Gram-Schmidt: Gedrag bij foutafvlakking
% Zie verslag


%% DEELVRAAG 7 Gram-Schmidt: Orthogonaliteit Maat en grafieken
% We gebruiken nu de Q_kl en Q_mod matrices die zijn berekend in DEELVRAAG 2

% Parameters
m = 200;
n = 40;

% Matrix A genereren: A = Q * D * V
Q = orth(randn(m,n));         % Orthogonale matrix (m x n)
D = diag(2.^-(0:n-1));        % Diagonale matrix met entries 2^-(i-1)
V = eye(n); V(1,:) = 1;       % V: eerste rij vol enen, rest identiek
A = Q * D * V;                % A = QDV

% Klassieke en Gewijzigde Gram-Schmidt toepassen
[Q_kl, ~] = klGS(A);
[Q_mod, ~] = modGS(A);

% Referentie-Q om fout te berekenen
Q_orig = Q;

% Foutnormen ||q̃j - qj||_2
error_kl = vecnorm(Q_kl - Q_orig, 2, 1);
error_mod = vecnorm(Q_mod - Q_orig, 2, 1);

% Orthogonaliteitsmaat: ||I - Q̃_j^T Q̃_j||_2 voor j = 1..n
orth_kl = arrayfun(@(j) norm(eye(j) - Q_kl(:,1:j)' * Q_kl(:,1:j), 2), 1:n);
orth_mod = arrayfun(@(j) norm(eye(j) - Q_mod(:,1:j)' * Q_mod(:,1:j), 2), 1:n);

% Plotten
figure;
semilogy(1:n, error_kl, '.-', 'DisplayName', 'Fout Klassieke GS', 'LineWidth', 1.5);
hold on;
semilogy(1:n, error_mod, '.-', 'DisplayName', 'Fout Gewijzigde GS', 'LineWidth', 1.5);
semilogy(1:n, orth_kl, '--', 'DisplayName', 'Orthogonaliteit Klassieke GS', 'LineWidth', 1.5);
semilogy(1:n, orth_mod, '--', 'DisplayName', 'Orthogonaliteit Gewijzigde GS', 'LineWidth', 1.5);

xlabel('Kolomindex $j$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Waarde (logaritmische schaal)', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Vergelijking van Fout en Orthogonaliteit}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 10, 'Interpreter', 'latex');
grid on;


%% DEELVRAAG 8 Gram-Schmidt: Genereren van matrices en conditiegetal berekenen
% Parameters
m = 200;
n = 50;
k_waardes = (0:16)'; 
conditie_getallen = zeros(size(k_waardes));

for i = 1:length(k_waardes)
    k = k_waardes(i);

    % Genereer een random orthogonale matrix Q
    Q = orth(randn(m, n));
    
    % Genereer een orthogonale vierkante matrix V
    V = orth(randn(n, n));
    
    % Genereer de diagonaalmatrix D
    D = diag(2 .^ linspace(0, k, n));
    
    % Matrix A = Q * D * V
    A = Q * D * V;
    
    % Bereken het conditiegetal
    conditie_getallen(i) = cond(A);
end

resultaat_tabel = table(k_waardes, conditie_getallen, 'VariableNames', {'k', 'Conditiegetal'});

disp(resultaat_tabel);

%% DEELVRAAG 9 Gram-Schmidt: Fout op de orthogonaliteit plotten voor drie algoritmen

% Parameters
m = 200;
n = 50;
k_max = 50;

k_values = 1:k_max;  % bepaalt hoe scherp de singular values afnemen
error_klassieke = zeros(size(k_values));
error_gewijzigde = zeros(size(k_values));
error_herhaalde = zeros(size(k_values));
kappa_values = zeros(size(k_values));

% Machineprecisie
epsmach = eps;

for k = k_values
    % Genereer vaste orthogonale matrices
    Q_base = orth(randn(m, n));
    V_base = orth(randn(n, n));

    % Maak D afhankelijk van k maar hou afmetingen gelijk
    D = diag(2 .^ linspace(0, -k, n));  % afname afhankelijk van k
    A = Q_base * D * V_base;  % volledige matrix A

    % Volledige orthogonalisatie
    [Q_kl, ~] = klGS(A);
    [Q_mod, ~] = modGS(A);
    [Q_her, ~] = herGS(A);

    % Bereken orthogonaliteitsfout op alle kolommen
    error_klassieke(k) = norm(eye(n) - Q_kl' * Q_kl, 'fro');
    error_gewijzigde(k) = norm(eye(n) - Q_mod' * Q_mod, 'fro');
    error_herhaalde(k) = norm(eye(n) - Q_her' * Q_her, 'fro');

    % Bepaal cond(A)
    kappa_values(k) = cond(A);
end

% Referentiecurves
ref_r1 = epsmach * kappa_values;
ref_r2 = epsmach * kappa_values.^2;
ref_r3 = epsmach * kappa_values.^0;

% Plot
figure;
semilogy(k_values, error_klassieke, '.-', 'LineWidth', 1.5, 'DisplayName', 'Klassieke GS');
hold on;
semilogy(k_values, error_gewijzigde, '.-', 'LineWidth', 1.5, 'DisplayName', 'Gewijzigde GS');
semilogy(k_values, error_herhaalde, '.-', 'LineWidth', 1.5, 'DisplayName', 'Herhaalde GS');
semilogy(k_values, ref_r1, '--k', 'LineWidth', 1.5, 'DisplayName', '$\epsilon_{\mathrm{mach}} \cdot \kappa(A)$');
semilogy(k_values, ref_r2, ':k', 'LineWidth', 1.5, 'DisplayName', '$\epsilon_{\mathrm{mach}} \cdot \kappa(A)^2$');
semilogy(k_values, ref_r3, '-.k', 'LineWidth', 1.5, 'DisplayName', '$\epsilon_{\mathrm{mach}} \cdot \kappa(A)^0$');
xlabel('Waarde $k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Fout op orthogonaliteit $\|I - \widetilde{Q}^\top \widetilde{Q}\|$', 'Interpreter', 'latex', 'FontSize', 12);
title('\textbf{Fout op orthogonaliteit bij volledige matrix $A$}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 10, 'Interpreter', 'latex');
grid on;


%% DEELVRAAG 10 Gram-Schmidt: adaptieve herhaalde GS
% Zie verslag


%% DEELVRAAG 11 Gram-Schmidt: Verifiëren van de orthogonaliteit van Q op basis van de gegeven matrix Amat

% Laad de matrix Amat
load('Amat.mat');  % Dit laadt de matrix Amat uit het bestand

% Voer de drie methoden uit (Klassieke, Gewijzigde, en Herhaalde GS)
[Q_kl, ~] = klGS(Amat);
[Q_mod, ~] = modGS(Amat);
[Q_herhaalde, ~] = herGS(Amat);

% Berekenen fout op de orthogonaliteit voor elk van de drie methoden
E_kl = Q_kl' * Q_kl - eye(size(Q_kl, 2));  % Fout voor Klassieke GS
E_mod = Q_mod' * Q_mod - eye(size(Q_mod, 2));  % Fout voor Gewijzigde GS
E_herhaalde = Q_herhaalde' * Q_herhaalde - eye(size(Q_herhaalde, 2));  % Fout voor Herhaalde GS

% Gebruik logaritme om de grootte-orde van de fouten te bepalen
log_E_kl = log10(abs(E_kl));
log_E_mod = log10(abs(E_mod));
log_E_herhaalde = log10(abs(E_herhaalde));

% Plot de heatmaps
figure('Position', [100, 100, 1000, 400]);  % Maak de figuur groter (breder)

subplot(1,3,1);
heatmap(log_E_kl);
title('Klassieke GS Fout');
colorbar; % Voeg kleurenschaal toe

subplot(1,3,2);
heatmap(log_E_mod);
title('Gewijzigde GS Fout');
colorbar; % Voeg kleurenschaal toe

subplot(1,3,3);
heatmap(log_E_herhaalde);
title('Herhaalde GS Fout');
colorbar; % Voeg kleurenschaal toe


%% DEELVRAAG 1 Splines: KKB met een kubische spline 
% zie kkb_cubespline.m , bsplineBasis.m en deBoor.m


%% DEELVRAAG 2 Splines: Correctheid algoritme van Boor

% Testdata
x = linspace(0, 10, 20)';     % meetpunten (r = 20)
b = sin(x);                  % functiewaarden

% Kies een knopenrij. Zorg dat de lengte voldoet aan: length(t) = nb + k + 1.
% Voor kubische splines (k = 3) kiezen we bijv.  nb = 10 basisfuncties.
k = 3;                       % graad (kubisch)
nb = 10;
m = nb + k + 1;

% Bepaal de knopenrij t zodat er k elementen voor min(x) en  k elementen na
% max(x) zijn.
h = (max(x) - min(x)) / (nb - k); 
t_start = min(x) - (k)*h;
t_end = max(x) + (k)*h ;
t = linspace(t_start, t_end, m);  % lineair verdeelde knopen met extra marge

% Bepaal evaluatiepunten en evalueer de spline met de Boor-algoritme
x_eval = linspace(min(x), max(x), 1000);

% Gebruik de functie kkb_cubespline voor de berekeningen
z = kkb_cubespline(t, x, b, x_eval);

% Plot de kleinstekwadratenbenadering en de meetpunten
figure;
plot(x_eval, z, 'b-', 'LineWidth', 1.5); hold on;
plot(x, b, 'ro', 'MarkerFaceColor', 'r');
plot(x_eval, sin(x_eval), 'k--', 'LineWidth', 2);  % Originele functie als gestreepte lijn
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Kleinstekwadratenbenadering met kubische B-splines}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Spline benadering', 'Meetpunten', 'Originele functie', 'Interpreter', 'latex', 'FontSize', 10);
grid on;

% Als illustratie: plot enkele individuele B-splines.
figure;
hold on;
% Zet coëfficiënten apart; voor elke basisfunctie activeren we telkens één coëfficiënt
for i = 1:nb
    c_single = zeros(nb, 1);
    c_single(i) = 1;
    z_single = zeros(size(x_eval));
    for j = 1:length(x_eval)
        z_single(j) = deBoor(x_eval(j), t, c_single, k);
    end
    plot(x_eval, z_single, 'LineWidth', 2);
end
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
 ylabel('$N_{i, k+1}(x)$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Enkele individuele B-splines (kubisch)}', 'Interpreter', 'latex', 'FontSize', 14);
legend(arrayfun(@(i) sprintf('Basis %d', i), 1:nb, 'UniformOutput', false), 'Interpreter', 'latex', 'FontSize', 10);
grid on;


%% DEELVRAAG 3 Splines: Belang van knopenligging bij B-spline benadering
% --- Meetpunten en functie ---
x = linspace(1, 9, 50)';   % Gebruik 50 meetpunten voor een gedetailleerdere benadering
% Definieer een "complexe" sinc-functie met lagere frequentie:
func_benader = @(y) (sin(2 * pi * (y - 5))) ./ (pi * (y - 5));
b = func_benader(x); % Functiewaarden

k = 3;      % Kubische B-splines
nb = 20;    % Meer basisfuncties voor een nauwkeurigere benadering

% 1. Gelijk verdeelde knopen
t1_internal = linspace(min(x), max(x), nb - k + 1);
t1 = [repmat(t1_internal(1), 1, k), t1_internal, repmat(t1_internal(end), 1, k)];
x_eval1 = linspace(t1(k+1), t1(end-k), 1000);
z1 = kkb_cubespline(t1, x, b, x_eval1);

% 2. Knopen geclusterd in het midden (middenclustering)
% We willen de interne knopen concentreren in het midden (ongeveer rond x=5).
% Definieer een uniforme parameter v in [0,1]:
v = linspace(0, 1, nb - k + 1);
% Pas een transformatie toe die de incrementen in het midden verkleint.
d = 0.1; % Een kleine parameter; grotere d geeft sterkere clustering, maar behoud d zodanig dat de functie monotoon blijft.
u = v - d * sin(2*pi*(v-0.5));
% Normaliseer u naar het interval [0,1]:
u = (u - min(u)) / (max(u) - min(u));
% Schaal naar het domein van x:
t2_internal = u * (max(x) - min(x)) + min(x);
t2 = [repmat(t2_internal(1), 1, k), t2_internal, repmat(t2_internal(end), 1, k)];
x_eval2 = linspace(t2(k+1), t2(end-k), 1000);
z2 = kkb_cubespline(t2, x, b, x_eval2);

% Visualisatie
figure;

% Subplot 1: Spline met gelijk verdeelde knopen
subplot(2,1,1);
plot(x_eval1, z1, 'b-', 'LineWidth', 1.5); hold on;
plot(x, b, 'ro', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
% Originele functie: bereken met de formule:
plot(x_eval1, func_benader(x_eval1), 'k--', 'LineWidth', 2);
ylim([-1.5, 3]);
title('\textbf{Spline met gelijk verdeelde knopen}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Spline benadering', 'Meetpunten', 'Originele functie', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);  
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;

% Subplot 2: Spline met knopen geclusterd in het midden
subplot(2,1,2);
plot(x_eval2, z2, 'm-', 'LineWidth', 1.5); hold on;
plot(x, b, 'ro', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
plot(x_eval2, func_benader(x_eval2), 'k--', 'LineWidth', 2);
ylim([-1.5, 3]);
title('\textbf{Spline met knopen geclusterd in het midden}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Spline benadering', 'Meetpunten', 'Originele functie', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);  
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;

fout_1 = abs(z1 - func_benader(x_eval1));
fout_2 = abs(z2 - func_benader(x_eval2));
figure;
plot(x_eval1, fout_1, 'b-', 'LineWidth', 1.5); hold on;
plot(x_eval2, fout_2, 'm-', 'LineWidth', 1.5);
title('\textbf{Foutnormen van de benaderingen}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Foutnorm', 'Interpreter', 'latex', 'FontSize', 13);
legend('Gelijke verdeling', 'Clustering in het midden', 'Interpreter', 'latex', 'FontSize', 10);
grid on;

% Knooppuntvisualisatie
figure;
hold on;
y_offset_1 = 2.5;
y_offset_2 = 2.6;
plot(t1, y_offset_1*ones(size(t1)), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(t2, y_offset_2*ones(size(t2)), 'mo', 'MarkerSize', 8, 'LineWidth', 1.5);
ylim([1, 4.1]);
yticks([]);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Ligging van de knopenrijen}', 'Interpreter', 'latex', 'FontSize', 14);
legend('Gelijke verdeling', 'Clustering in het midden', 'Interpreter', 'latex', 'FontSize', 10);
grid on;


%% DEELVRAAG 4 Splines: Effect van ruis op splinebenadering

% Doelfunctie (ruisvrij)
x = linspace(-1, 1, 200)';
b = sin(20*x) ./ (100*x.^2 + 5);

% Toevoegen van ruis
rng(0); % voor reproduceerbaarheid
b_ruis = b + 0.04 * randn(size(x));

% Instellingen
k = 3; % kubische B-spline
max_nb = 42; % maximum aantal basisfuncties
min_nb = 6; % minimum aantal basisfuncties
nb_range = min_nb:max_nb;

err_true = zeros(size(nb_range));
err_ruis = zeros(size(nb_range));

% Evaluatie over verschillende aantallen knopen
for i = 1:length(nb_range)
    nb = nb_range(i);
    % Lineair verdeelde knopen
    h = (max(x) - min(x)) / (nb - k); 
    t_start = min(x) - (k)*h;
    t_end = max(x) + (k)*h ;
    t = linspace(t_start, t_end, nb + k + 1); 
    
    % Spline op ruisdata
    z = kkb_cubespline(t, x, b_ruis, x);
    
    % Fouten
    err_true(i) = norm(b - z);
    err_ruis(i) = norm(b_ruis - z);
end

x_knopen = nb_range - k + 1;  % correcte maat voor het aantal (interne) knopen

% ==== Plot van foutnormen ====
figure;
plot(x_knopen, err_true, 'b-o', 'LineWidth', 1.5); hold on;
plot(x_knopen, err_ruis, 'r-s', 'LineWidth', 1.5);
xlabel('Aantal (interne) knopen', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('Foutnorm', 'Interpreter', 'latex', 'FontSize', 13);
legend('$\|b - z\|$ (ware fout)', '$\|b_{\mathrm{ruis}} - z\|$ (residuele fout)', ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
title('\textbf{Effect van knoopaantal op de benaderingskwaliteit}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% ==== Beste benaderingen ====
[~, idx_true] = min(err_true);
[~, idx_ruis] = min(err_ruis);

nb_true = nb_range(idx_true);
nb_ruis = nb_range(idx_ruis);

% Knopenrij voor true
h_true = (max(x) - min(x)) / (nb_true - k); 
t_start_true = min(x) - (k)*h_true;
t_end_true = max(x) + (k)*h_true ;
t_true = linspace(t_start_true, t_end_true, nb_true + k + 1); 

% Knopenrij voor ruis
h_ruis = (max(x) - min(x)) / (nb_ruis - k); 
t_start_ruis = min(x) - (k)*h_ruis;
t_end_ruis = max(x) + (k)*h_ruis ;
t_ruis = linspace(t_start_ruis, t_end_ruis, nb_ruis + k + 1); 

z_best_true = kkb_cubespline(t_true, x, b_ruis, x);
z_best_ruis = kkb_cubespline(t_ruis, x, b_ruis, x);

% Plot: Beste benadering voor ||f - z||
figure;
plot(x, b, 'k--', 'LineWidth', 1.5); hold on;
plot(x, z_best_true, 'b-', 'LineWidth', 1.5);
legend('Originele functie', 'Beste spline $(\|b - z\|)$', ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northwest'); % Legenda links boven
title(['\textbf{Beste benadering voor } $\mathbf{\|b - z\|}$ \textbf{(' num2str(nb_true) ' basisfuncties)}'], ...
      'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;

% Plot: Beste benadering voor ||f_{ruis} - z||
figure;
plot(x, b, 'k--', 'LineWidth', 1.5); hold on;
plot(x, z_best_ruis, 'r-', 'LineWidth', 1.5);
legend('Originele functie', 'Beste spline $(\|b_{\mathrm{ruis}} - z\|)$', ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northwest'); % Legenda links boven
title(['\textbf{Beste benadering voor } $\mathbf{\|b_{\mathrm{ruis}} - z\|}$ \textbf{(' num2str(nb_ruis) ' basisfuncties)}'], ...
      'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;