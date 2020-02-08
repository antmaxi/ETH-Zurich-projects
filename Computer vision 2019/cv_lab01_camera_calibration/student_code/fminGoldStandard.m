%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xy_normalized: 3xn
% XYZ_normalized: 4xn

function f = fminGoldStandard(pn, xy_normalized, XYZ_normalized)

%reassemble P

%compute reprojection error
s = size(xy_normalized);
n = s(1,2);
error = 0;


for i = 1:n
    XYZ_normalized([1,2,3],i) = XYZ_normalized([1,2,3],i)/XYZ_normalized(4,i);
    XYZ_normalized(4,i) = 1;
end
XYZ_P = pn*XYZ_normalized;
for i = 1:n
    error = error + norm(xy_normalized([1,2],i)/xy_normalized(3,i) - XYZ_P([1,2],i)/XYZ_P(3,i))^2;
end
%compute cost function value
f = sqrt(error/n);
end