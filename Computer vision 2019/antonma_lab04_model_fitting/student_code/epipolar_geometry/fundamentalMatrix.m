% Compute the fundamental matrix using the eight point algorithm
% Input
% 	x1s, x2s 	Point correspondences 3xN
%
% Output
% 	Fh 		Fundamental matrix with the det F = 0 constraint
% 	F 		Initial fundamental matrix obtained from the eight point algorithm
%
function [Fh, F] = fundamentalMatrix(x1s, x2s)
    % Normalize input points
    [nx1s, T1] = normalizePoints2d(x1s);
    [nx2s, T2] = normalizePoints2d(x2s);
    % Produce matrix of equation in order to obtain later fundamental
    % matrix
    u = nx1s(1, :)';
    v = nx1s(2, :)';
    u1 = nx2s(1, :)';
    v1 = nx2s(2, :)';
    A = [u.*u1, u.*v1, u, ...
         v.*u1, v.*v1, v, ...
           u1,   v1, ones(size(u,1),1)];
    [U,S,V] = svd(A);
    f = V(:, end);
    % Get fundamental matrix
    F = reshape(f, 3, 3);
    [U1,S1,V1] = svd(F);
    % Normalized fundamental matrix
    F = T2'*F*T1;
    % Singularize fundamental matrix
    S1(3,3) = 0;
    Fh = T2'*U1*S1*V1'*T1;
end