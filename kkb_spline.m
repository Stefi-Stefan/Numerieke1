function Z = kkb_spline(t, x, f, xplot)
% Deze functie maakt een kleinstekwadratenbenadering met B-splines
% en evalueert die benadering in de punten van xplot.
%
%   Z = kkb_spline(t, x, f, xplot)
%
%   Input:
%       t     - kolomvector met knopen(lengte n+2k+1) 
%       x     - kolomvector met meetpunten (lengte R)
%       f     - kolomvector met functiewaarden in x (lengte R)
%       xplot - rijvector met evaluatiepunten (lengte N)
%
%   Output:
%       Z     - rijvector met splinewaarden in xplot

% Zet alles in de juiste vorm.
t = t(:);
x = x(:);
f = f(:);
xplot = xplot(:).';

% In deze opdracht werken we met kubische splines.
k = 3;

% Aantal basisfuncties.
nb = length(t) - k - 1;

% Bouw de matrix M van het overgedetermineerde stelsel M*c = f.
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

    % Enkel k+1 basisfuncties zijn niet nul in xr.
    Nloc = basisFuns(span, xr, k, t);
    indices = (span-k):span;
    M(r, indices) = Nloc;
end

% Los het kleinstekwadratenprobleem op.
c = M \ f;

% Evalueer de spline in de gevraagde punten met deBoor.
Z = zeros(1, length(xplot));
for j = 1:length(xplot)
    Z(j) = deBoor(xplot(j), t, c, k);
end

end

function span = findSpan(nb, k, x, t)
% Bepaal de span-index zodat t(span) <= x < t(span+1).

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
% Bereken de lokale niet-nul basisfuncties in punt x.

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
