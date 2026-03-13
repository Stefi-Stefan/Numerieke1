% Gewijzijgd Gram-Schmidt Algoritme
% Zie cursustekst sectie 2.3.2 'Algoritme 2'
function [Q, R] = modGS(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    
    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * v;
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end