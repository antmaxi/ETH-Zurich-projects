%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xy: size 2xn
% XYZ: size 3xn 

function [P, K, R, t, error] = runDLT(xy, XYZ)
disp('data');
disp(xy);
disp('000');
disp(XYZ);
% normalize 
[xy_normalized, XYZ_normalized, T, U] = normalization(xy, XYZ)

%compute DLT with normalized coordinates
[Pn] = dlt(xy_normalized, XYZ_normalized);
%denormalize projection matrix
Pn = [Pn([1,2,3,4],1),Pn([5,6,7,8],1),Pn([9,10,11,12],1)].';
%disp('Pn');
%disp(Pn);
P = T^(-1)*Pn*U;
%factorize projection matrix into K, R and t
[K, R, t] = decompose(P);
%disp('P');
%disp(P);

%compute average reprojection error
error = 0;
s = size(xy);
n = s(1, 2);
%disp('xxx');
%disp(xy_normalized(:,1));
XYZ_hom = homogenization(XYZ);
XYZ_P = P*XYZ_hom;
for i = 1:n
    error = error + norm(xy(:,i) - XYZ_P([1,2],i)/XYZ_P(3,i))^2;
end
error = sqrt(error/n);
disp('error DLT');
disp(error);
end