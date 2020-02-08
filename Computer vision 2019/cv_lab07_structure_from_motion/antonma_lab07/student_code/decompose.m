%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%

function [K, R, t] = decompose(P)
%Decompose P into K, R and t using QR decomposition

% Compute R, K with QR decomposition such M=K*R 
[R, K] = qr((P([1,2,3], [1,2,3]))^(-1));
R = R^(-1);
K = K^(-1);
% Compute camera center C=(cx,cy,cz) such P*C=0 
[U,S,V] = svd(P);
disp(size(V));
C = V(:, 4);
C = C/C(4);
% normalize K such K(3,3)=1
R = R*K(3,3);
K = K/K(3,3);
% Adjust matrices R and Q so that the diagonal elements of K = intrinsic matrix are non-negative values and R = rotation matrix = orthogonal has det(R)=1
K = K*abs(det(R))^(1/3);
%K(2:)= -1*K(2:);
R = R/abs(det(R))^(1/3);
disp('K');
disp(K);
disp('R');
disp(R);
disp(det(R))


% Compute translation t=-R*C
t = -R*C([1,2,3]);
end