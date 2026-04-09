% Herhaalde Gram-Schmidt Algoritme
% Zie cursustekst sectie 2.3.2 bij paragraaf met 'Stapsgewijse variant' vooraan. (P 36)
function [Q, R] = herGS(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    S = zeros(m, n);
  
    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:,j);
            v = v - R(i, j) * Q(:, i);
        end
        w = v;
        for i = 1:j-1
            S(i, j) = Q(:, i)' * w;
            v = v - S(i, j) * Q(:, i);
            R(i, j) = R(i, j) + S(i, j);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end