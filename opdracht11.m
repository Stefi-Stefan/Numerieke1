% structuur van matrix M en bespreking van efficiente QR-decompositie
clear; clc; close all;

S = load('exercise10.mat');
x = S.x(:);
f = S.y(:);
t1 = S.t1(:);
t2 = S.t2(:);

% kies een van de twee knopenrijen, hier nemen we t2 (beste fit oef 10)
t = t2;
k = 3; % kubische splines

% bouw matrix M van stelsel M*c = f
M = bouwM(t, x, k);

% zet om naar sparse voor duidelijke spy-plot en efficiente QR
Ms = sparse(M);

% visualiseer de structuur van M
figure;
spy(Ms);
title('Sparsity-structuur van M (opgave 11)');
xlabel('kolomindex');
ylabel('rijindex');

[Rm, Cm] = size(Ms); % Rm = aantal rijen (meetpunten R)
% Cm = aantal kolommen (basisfuncties n+k)
nzM = nnz(Ms); % nnz = Number of Non-Zero elements
% telt alle elementen die niet 0 zijn
densiteit = nzM / (Rm * Cm);  % verhouding van niet-nul elementen ten opzichte van de hele matrix
rowNnz = full(sum(Ms ~= 0, 2));
maxNnzRow = max(rowNnz); % maximale aantal actieve basisfuncties in 1 meetpunt (theoretisch max k+1)
minNnzRow = min(rowNnz); %minimale aantal actieve basisfuncties

% benaderde bandbreedte: uiterste kolommen met niet-nul per rij
onderBand = 0;
bovenBand = 0;
for r = 1:Rm
    cols = find(Ms(r, :));% kolom-indices van alle niet-nul elementen op rij r
    if ~isempty(cols)
        onderBand = max(onderBand, r - cols(1));% afstand diagonaal tot het eerste niet-nul element
        bovenBand = max(bovenBand, cols(end) - r);% afstand diagonaal tot het laatste niet-nul element
    end
end

fprintf('\nSTRUCTUUR VAN M\n');
fprintf('afmetingen M: %d x %d\n', Rm, Cm);
fprintf('aantal niet-nul elementen: %d\n', nzM);
fprintf('densiteit: %.3f %%\n', 100*densiteit);
fprintf('min aantal niet-nul per rij: %d\n', minNnzRow);
fprintf('max aantal niet-nul per rij: %d\n', maxNnzRow);
fprintf('verwacht voor kubisch: maximaal k+1 = %d niet-nul per rij\n', k+1);
fprintf('geschatte onderband: %d\n', onderBand);
fprintf('geschatte bovenband: %d\n', bovenBand);

% QR op sparse matrix (economy)
[Q, R] = qr(Ms, 0);

figure;
spy(R);
title('Sparsity-structuur van R bij sparse QR van M');
xlabel('kolomindex');
ylabel('rijindex');

fprintf('\nEFFICIENTE QR-BEREKENING\n');
fprintf('1) M is schaars en bandvormig (door lokale steun van B-splines).\n');
fprintf('2) Gebruik sparse QR in plaats van dense QR : qr(sparse(M),0).\n');
fprintf('3) Werk met transformaties op f (Q''*f) zonder expliciet volledige Q op te slaan.\n');
fprintf('4) Zo beperk je rekentijd en geheugen, zeker voor grote datasets.\n');

% kleinstekwadratenbenadering oplossen via QR en vergelijken met backslash
c_qr = R \ (Q' * f);
c_bs = Ms \ f;
relVerschil = norm(c_qr - c_bs, 2) / (norm(c_bs, 2) + eps);
fprintf('relatief verschil tussen QR-oplossing en backslash: %.3e\n', relVerschil);

function M = bouwM(t, x, k)
% bouw matrix M met M(r,i) = N_i,k+1(x_r)

nb = length(t) - k - 1;
R = length(x);
M = zeros(R, nb);

for r = 1:R
    xr = x(r);

    if xr < t(k+1) || xr > t(nb+1)
        continue;
    end

    if xr == t(nb+1)
        span = nb;
    else
        span = findSpan(nb, k, xr, t);
    end

    Nloc = basisFuns(span, xr, k, t);
    idx = (span-k):span;
    M(r, idx) = Nloc;
end
end

function span = findSpan(nb, k, x, t)
% bepaal span-index zodat t(span) <= x < t(span+1)

low = k + 1;
high = nb + 1;
mid = floor((low + high) / 2);
while ~(x >= t(mid) && x < t(mid+1))
    if x < t(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end
span = mid;
end

function Nloc = basisFuns(i, x, k, t)
% bereken lokale niet-nul basisfuncties in punt

Nloc = zeros(1, k+1);
Nloc(1) = 1;
left = zeros(1, k);
right = zeros(1, k);

for j = 1:k
    left(j) = x - t(i+1-j);
    right(j) = t(i+j) - x;
    saved = 0;

    for r = 1:j
        temp = Nloc(r) / (right(r) + left(j+1-r));
        Nloc(r) = saved + right(r) * temp;
        saved = left(j+1-r) * temp;
    end

    Nloc(j+1) = saved;
end
end
