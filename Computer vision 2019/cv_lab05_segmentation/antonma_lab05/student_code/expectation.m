% Input
% Image X - Lx3
% alpha Kx1
% mu - 3xK
% var - 3x3xK
%%%%%%%%%%%%%%%%
% Output
% P - LxK
function P = expectation(mu,var,alpha,X)

K = length(alpha);
% number of pixels
L = size(X,1);
P = zeros(L, K);
% prepare array of inversed covariance matrices and determinants
inverse = zeros(3,3,K);
sqr_det = zeros(K,1);
for i = 1:K
    inverse(:,:,i) = inv(var(:,:,i));
    sqr_det(i) = det(var(:,:,i))^0.5;
end
% Iterate over pixels
for i = 1:L
    % Normalization constant
    Z = 0;
    for j = 1:K
    % pdf of multivariate Gaussian
        P(i,j) = alpha(j)*(2*pi)^(-3/2)/sqr_det(j)*exp(-0.5*(X(i,:)'-mu(:,j))'*inverse(:,:,j)*(X(i,:)'-mu(:,j)));
        Z = Z + P(i,j);
    end
    P(i, :) = P(i, :)/Z;
end
end