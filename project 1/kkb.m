function C = kkb(x, y, F, n, m)
% Deze functie lost het 2D-kleinstekwadratenprobleem op met Legendrebasis.
%
%   C = kkb(x, y, F, n, m)
%
%   Input:
%       x - rijvector met x-waarden (lengte N)
%       y - rijvector met y-waarden (lengte M)
%       F - datamatrix van grootte M x N
%       n - maximale graad in x
%       m - maximale graad in y
%
%   Output:
%       C - coefficientenmatrix van grootte (m+1) x (n+1)
%           met C = (B^+) F (A^+)^T

x = x(:);
y = y(:);
A = get_leg_mtx(x, n); % dimensie N x (n+1)
B = get_leg_mtx(y, m); % dimensie M x (m+1)

%B * X = F dus X = B \ F geeft X = (B^+) * F
X = B \ F;

%A * Y = X' dus Y = A \ X' geeft Y = (A^+) * X'
Y = A \ X';
C = Y';

end
