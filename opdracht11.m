% Opgave 11
% Structuur van matrix M en bespreking van efficiente QR-decompositie.
clear; clc; close all;
tic

% Laad de data van oefening 10.
S = load('exercise10.mat');
x = S.x(:);
y = S.y(:);
t1 = S.t1(:);
t2 = S.t2(:);

% Kies een van de twee knopenrijen.
% Hier nemen we t2 (beste fit oef 10)
t = t2;

% Kubische splines.
k = 3;

% Bouw matrix M van stelsel M*c = y.
M = bouwM(t, x, k);

% Zet om naar sparse voor duidelijke spy-plot en efficiente QR.
Ms = sparse(M);

% Visualiseer de structuur van M.
figure;
spy(Ms);
title('Sparsity-structuur van M (opgave 11)');
xlabel('kolomindex');
ylabel('rijindex');



function M = bouwM(t, x, k)
% Bouw matrix M met M(r,i) = N_i,k+1(x_r).

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
% Bepaal span-index zodat t(span) <= x < t(span+1).

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
% Bereken lokale niet-nul basisfuncties in punt x.

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
toc