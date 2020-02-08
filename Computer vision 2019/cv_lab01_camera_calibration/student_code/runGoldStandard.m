%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xy: size 2xn
% XYZ: size 3xn 

function [P, K, R, t, error] = runGoldStandard(xy, XYZ)

%normalize data points
[xy_normalized, XYZ_normalized, T, U] = normalization(xy, XYZ);
xy_hom = homogenization(xy);
XYZ_hom = homogenization(XYZ);
%compute DLT with normalized coordinates
[pn] = dlt(xy_normalized, XYZ_normalized);
p1 = [pn([1,2,3,4],1),pn([5,6,7,8],1),pn([9,10,11,12],1)].';
P2 = T^(-1)*p1*U;

IMG_NAME = 'image001.jpg';
visualization_reprojected_points(xy, XYZ, P2, IMG_NAME);
%minimize geometric error to refine P_normalized
options = optimset('MaxFunEvals',50000,'MaxIter', 20000);
for i=1:20
    [p1] = fminsearch(@fminGoldStandard, p1, [options], xy_normalized, XYZ_normalized);
end

%denormalize projection matrix
P = T^(-1)*p1*U;
%factorize prokection matrix into K, R and t
[K, R, t] = decompose(P);

%compute average reprojection error
error = 0;
s = size(xy);
n = s(1, 2);
XYZ_hom = homogenization(XYZ);
XYZ_P = P*XYZ_hom;
for i = 1:n
    error = error + norm(xy(:,i) - XYZ_P([1,2],i)/XYZ_P(3,i))^2;
end
%compute cost function value
f = sqrt(error/n);
disp('GoldError')
disp(f);
end