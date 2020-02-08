% Input
% Image X - Lx3
% Probabilities P - LxK
%%%%%%%%%%%%%%%%
% Output
% mu - 3xK
% var - 3x3xK
% alpha Kx1
function [mu, var, alpha] = maximization(P, X)

K = size(P,2);
L = size(X,1);
alpha = zeros(K,1);
mu = zeros(3,K);
var = zeros(3,3,K);
for i = 1:K
    % Normalizing constant
    Z = 0;
    for j = 1:L
        Z = Z + P(j, i);
        mu(:, i) = mu(:, i) + X(j,:)'*P(j, i);
    end
    alpha(i) = Z/L;
    mu(:, i) = mu(:, i)/Z;
    for j = 1:L
        var(:,:,i) = var(:,:,i) + P(j, i)*(X(j,:)' - mu(:, i))*(X(j,:)' - mu(:, i))';
    end
    % normalize and add small epsilon to covariance matrix in order to 
    % preserve positive semi-definite property, which could break due to rounding
    var(:,:,i) = var(:,:,i)/Z + 10^(-10)*eye(3);
end