clear; clc; close all;

k = 3;% kubische splines dus k = 3.
t = [0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1];% de knopenrij is een open geklemde knopenrij met
% inwendige knopen 0.2, 0.4, 0.6 en 0.8.

nb = length(t) - k - 1; % length(t)=n+2k+1 aantal knooppunten, en
% het aantal basisfuncties is n+k = n+2k+1-k-1
xplot = linspace(0, 1, 800);
B = zeros(nb, length(xplot));

% bereken de basisfuncties door telkens 1 coefficient gelijk aan 1 te
% nemen, en de rest 0
for i = 1:nb
    c = zeros(nb, 1);
    c(i) = 1;
    B(i,:) = deboor(k, t, c, xplot);
end

figure;
plot(xplot, B, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('N_{i,k+1}(x)');
title(sprintf('B-splinebasis op [0,1], graad k = %d', k));
legend('N_{-3,4}','N_{-2,4}','N_{-1,4}','N_{0,4}','N_{1,4}','N_{2,4}','N_{3,4}','N_{4,4}', 'Location', 'eastoutside');

disp('Gekozen parameters voor oefening 8:');
disp(['k = ' num2str(k)]);
disp(['aantal basisfuncties = ' num2str(nb)]);
disp('knopenrij t =');
disp(t);

% controle: de som van de basisfuncties moet ongeveer 1 zijn
partitionError = max(abs(sum(B, 1) - 1));
disp(['maximale fout op partition of unity = ' num2str(partitionError)]);
