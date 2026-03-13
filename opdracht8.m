% Oefening 8
% Visualisatie van de volledige B-splinebasis op [0,1]
clear; clc; close all;

% Kubische splines, dus k = 3.
% De knopenrij is een open geklemde knopenrij met
% inwendige knopen 0.2, 0.4, 0.6 en 0.8.
k = 3;
t = [0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1];

nb = length(t) - k - 1;
xplot = linspace(0, 1, 800);
B = zeros(nb, length(xplot));

% Bereken de basisfuncties door telkens 1 coefficient gelijk aan 1 te nemen.
for i = 1:nb
    c = zeros(nb, 1);
    c(i) = 1;
    for j = 1:length(xplot)
        B(i,j) = deBoor(xplot(j), t, c, k);
    end
end

figure;
plot(xplot, B', 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('N_{i,k+1}(x)');
title(sprintf('B-splinebasis op [0,1], graad k = %d', k));
legend('N_{1,4}','N_{2,4}','N_{3,4}','N_{4,4}','N_{5,4}','N_{6,4}','N_{7,4}','N_{8,4}', 'Location', 'eastoutside');

disp('Gekozen parameters voor oefening 8:');
disp(['k = ' num2str(k)]);
disp(['aantal basisfuncties = ' num2str(nb)]);
disp('knopenrij t =');
disp(t);

% Controle: de som van de basisfuncties moet ongeveer 1 zijn.
partitionError = max(abs(sum(B, 1) - 1));
disp(['maximale fout op partition of unity = ' num2str(partitionError)]);
