function N = bsplineBasis(t, k, x)
% Deze functie berekent de globale B-spline basisfunctiewaarden bij een scalar x.
%
%   N = bsplineBasis(t, k, x) geeft een rijvector terug van lengte nb = length(t)-k-1,
%   waarbij elk element de waarde van de corresponderende B-spline basisfunctie van orde k+1 
%   is, berekend met de Cox-deBoor-recusie.
%   Zie cursustekst sectie 5.3.3 eigenschap 5.3
%
%   Input:
%       t - Knooppuntenrij
%       k - Graad van de B-spline
%       x - Evaluatiepunt
%
%   Output:
%       N - Rijvector van lengte nb met de waarde van de i-de B-spline basisfunctie.

% Controleer of x binnen het geldige interval ligt:
nb = length(t) - k - 1;
if x < t(k+1) || x > t(nb+1)
    N = zeros(1, nb);
    return;
end

% Speciale behandeling: als x exact gelijk is aan het laatste knooppunt, dan krijgt de
% laatste basisfunctie waarde 1.
if x == t(nb+1)
    N = zeros(1, nb);
    N(end) = 1;
    return;
end

% Bepaal de "span" index via een binaire zoekmethode
span = findSpan(nb, k, x, t);

% Bereken de lokale basisfuncties (Cox-deBoor)
N_loc = basisFuns(span, x, k, t);

% Vul de globale basisfunctievector in (de niet-nul waarden staan in de indices span-k+1 tot span+1)
N = zeros(1, nb);
globalIndices = (span - k):(span);
N(globalIndices) = N_loc;

end

% --- Helperfunctie: findSpan ---
function span = findSpan(nb, k, x, t)
% Bepaal de span-index zodat t(span) <= x < t(span+1).
% Neem aan dat x < t(nb+1) (de bovengrens) door de eerdere controle.
low = k+1; 
high = nb+1; % Dit komt overeen met t(nb+1)
mid = floor((low + high)/2);
while ~(x >= t(mid) && x < t(mid+1))
    if x < t(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high)/2);
end
span = mid;
end

% --- Helperfunctie: basisFuns ---
function N_loc = basisFuns(i, x, k, t)
% basisFuns Bereken de lokale niet-nul waarden van de B-spline basisfuncties.
% Geeft een vector N_loc van lengte k+1, zodat N_loc(j+1) overeenkomt met N_{i-k+j}(x)
N_loc = zeros(1, k+1);
N_loc(1) = 1;
left = zeros(1, k);
right = zeros(1, k);
for j = 1:k
    left(j) = x - t(i+1-j);
    right(j) = t(i+j) - x;
    saved = 0;
    for r = 1:j
        temp = N_loc(r) / (right(r) + left(j+1-r));
        N_loc(r) = saved + right(r)*temp;
        saved = left(j+1-r)*temp;
    end
    N_loc(j+1) = saved;
end
end
