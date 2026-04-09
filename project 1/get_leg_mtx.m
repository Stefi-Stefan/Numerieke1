function A = get_leg_mtx(x, n)
% Deze functie maakt de vandermondematrix met Legendreveeltermen.
%
%   A = get_leg_mtx(x, n)
%
%   Input:
%       x - vector met evaluatiepunten (lengte N)
%       n - maximale graad
%
%   Output:
%       A - matrix van grootte N x (n+1)
%           met A(i,k+1) = P_k(x_i)

x = x(:);
N = length(x);
A = zeros(N, n+1);

% P_0(x) = 1
A(:,1) = 1;
if n >= 1
    % P_1(x) = x
    A(:,2) = x;
end

%recursie eig: (k+1)P_{k+1}(x) = (2k+1)xP_k(x) - kP_{k-1}(x)
for k = 1:n-1
    A(:,k+2) = ((2*k+1)*x.*A(:,k+1) - k*A(:,k)) / (k+1);
end

end
