function s = deBoor(x, t, c, k)
% Deze functie berekent de splinewaarde s(x) door het algoritme van de Boor.
% Zie cursustekst sectie 5.3.4.1
%
%   s = deBoor(x, t, c, k) evalueert de spline
%       s(x) = sum_{i=1}^{nb} c(i) N_i,k+1(x)
%   met de B-spline basisfuncties N_i,k+1(x).
%
%   Input:
%       x - Evaluatiepunt (scalar)
%       t - Knooppuntenrij (rij- of kolomvector)
%       c - Coefficienten van de spline (rij- of kolomvector)
%       k - Graad van de spline
%
%   Output:
%       s - Splinewaarde in x.

% Bepaal het aantal basisfuncties:
nb = length(t) - k - 1;

% Bepaal de juiste span-index i in de knopenrij t waarin x ligt.
if x == t(nb+1)
    i = nb;
else
    i = findSpanForDeBoor(nb, k, x, t);
end

% Zet de eerste set d's op basis van de coefficienten:
d_vec = zeros(1, k+1);
for j = 0:k
    % d_vec(j+1) correspondeert met c(i - k + j)
    d_vec(j+1) = c(i - k + j);
end

% Recursieve update van d_vec (algoritme van de Boor)
for r = 1:k
    for j = k:-1:r
        alpha = (x - t(i - k + j)) / (t(j + i - r + 1) - t(i - k + j));
        d_vec(j+1) = (1 - alpha) * d_vec(j) + alpha * d_vec(j+1);
    end
end

% De uiteindelijke splinewaarde:
s = d_vec(k+1);

end

% --- Helperfunctie: findSpanForDeBoor ---
function span = findSpanForDeBoor(nb, k, x, t)
% Bepaal de span-index i voor de deBoor-algoritme zodanig dat:
%   t(i) <= x < t(i+1)

% Speciaal geval: als x exact gelijk is aan de laatste knoop,
% return laatste geldige span
if abs(x - t(end)) < 1e-12
    span = nb;
    return;
end

low = k+1;
high = nb+1;
mid = floor((low+high)/2);
while ~(x >= t(mid) && x < t(mid+1))
    if x < t(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
span = mid;
end
