%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xyn: 3xn
% XYZn: 4xn

function [P_normalized] = dlt(xyn, XYZn)
%computes DLT, xy and XYZ should be normalized before calling this function

% 1. For each correspondence xi <-> Xi, computes matrix Ai
for i = 1:6
    V = XYZn(:,i);
    v = xyn(:,i);
    if i == 1
        A = [-V.',0,0,0,0,v(1)*V.';
            0,0,0,0,-V.',v(2)*V.'];
    else
        A = [A;
            -V.',0,0,0,0,v(1)*V.';
            0,0,0,0,-V.',v(2)*V.'];
    end
end
         
% 2. Compute the Singular Value Decomposition of A
[U,S,V] = svd(A);
% 3. Compute P_normalized (=last column of V if D = matrix with positive
% diagonal entries arranged in descending order)
P_normalized = V(:, 12);
end
