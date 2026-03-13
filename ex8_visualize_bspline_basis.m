% Oefening 8 - Visualiseer de volledige B-splinebasis op [0, 1]
% Nieuwe implementatie voor deze opdracht (zonder hergebruik van oude scripts).
clear; clc; close all;

% Gekozen opstelling (expliciet vermeld zoals gevraagd):
% - Kubische splines: k = 3
% - Open geklemde knopenrij op [0,1]
% - Interne knopen uniform op 0.2, 0.4, 0.6, 0.8
k = 3;
t = [0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1];

nb = length(t) - k - 1;    % aantal basisfuncties
xplot = linspace(0, 1, 800);
B = zeros(nb, length(xplot));

% Bouw elke basisfunctie op met een eenheidsvector als coefficienten.
for i = 1:nb
    c = zeros(nb, 1);
    c(i) = 1;
    for j = 1:length(xplot)
        B(i, j) = deBoor_ex8(xplot(j), t, c, k);
    end
end

figure('Color', 'w');
plot(xplot, B', 'LineWidth', 1.4);
grid on;
xlabel('x');
ylabel('N_{i,k+1}(x)');
title(sprintf('B-splinebasis op [0,1], graad k = %d', k));
legend(arrayfun(@(i) sprintf('N_{%d,%d}', i, k+1), 1:nb, 'UniformOutput', false), ...
    'Location', 'eastoutside');

fprintf('Instellingen oefening 8:\n');
fprintf('  Graad k = %d\n', k);
fprintf('  Knopenrij t = [ ');
fprintf('%.1f ', t);
fprintf(']\n');
fprintf('  Aantal basisfuncties = %d\n', nb);

% Optionele controles op basiseigenschappen op [0,1]:
partitionError = max(abs(sum(B, 1) - 1));
nonNegativityMin = min(B(:));
fprintf('  Maximale fout op partition of unity: %.3e\n', partitionError);
fprintf('  Minimale basiswaarde (moet >= 0 zijn): %.3e\n', nonNegativityMin);
