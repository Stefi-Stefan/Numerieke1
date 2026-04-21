function Phi = build_poly_features(x1, x2, n)
%BUILD_POLY_FEATURES Bouwt monomiale features voor model (3).
%   Phi = [1, x1, x2, x1^2, x2^2, ..., x1^n, x2^n]

N = numel(x1);
Phi = ones(N, 2*n + 1);

for k = 1:n
    Phi(:, 2*k)   = x1.^k;
    Phi(:, 2*k+1) = x2.^k;
end

end
