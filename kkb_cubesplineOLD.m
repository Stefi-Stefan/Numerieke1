function z = kkb_cubespline(t, x, b, xeval)
% Deze functie berekent de kleinstekwadratenbenadering met kubische B-splinefuncties.
%   Input:
%       t     - Knooppuntenrij
%       x     - Abscissen-vector met lengte r.
%       b     - Functiewaarden bij de punten x. kolomvector van lengte r.
%       xeval - Evaluatiepunten-vector met lengte N
%
%   Output:
%       z     - vector van lengte N met functiewaarden van de splinebenadering in de punten xecal

    k = 3;                      % Graad (kubisch)
    nb = length(t) - k - 1;     % Aantal B-spline basisfuncties
    r = length(x);            % Aantal meetpunten
    
    % Bouw de designmatrix M (r x nb)
    M = zeros(r, nb);
    for j = 1:r
        M(j, :) = bsplineBasis(t, k, x(j));
    end
    
    % Bereken de kleinstekwadratenoplossing voor de coëfficiënten c (overgedetermineerd stelsel)
    c = M \ b;
    
    % Evalueer de spline bij de punten in xeval met behulp van de Boor-algoritme
    z = zeros(size(xeval));
    for j = 1:length(xeval)
        z(j) = deBoor(xeval(j), t, c, k);
    end
end
