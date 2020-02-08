%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xy: size 2xn
% XYZ: size 3xn 

function [P, K, R, t, error] = runDLT(xy, XYZ)
disp('data');
disp(xy);
disp('000');
disp(XYZ);
% normalize 
% x3_sel_hom - 3xN
xyh = makehomogeneous(xy);
XYZh = makehomogeneous(XYZ);
[xy_normalized, T] = normalise2dpts(xyh);
[XYZ_normalized, U] = normalise3dpts(XYZh);
xyinh = makeinhomogeneous(xy_normalized);
XYZinh = makeinhomogeneous(XYZ_normalized);
%compute DLT with normalized coordinates
[Pn] = dlt(xy_normalized, XYZ_normalized);
%disp('Pn');
%disp(Pn);
%denormalize projection matrix
Pn = [Pn([1,2,3,4],1),Pn([5,6,7,8],1),Pn([9,10,11,12],1)].';
%disp('Pn');
%disp(Pn);
P = T^(-1)*Pn*U;
%factorize projection matrix into K, R and t
[K, R, t] = decompose(P);
disp('P');
disp(P);

%compute average reprojection error
error = 0;
s = size(xy);
n = s(1, 2);
%disp('xxx');
%disp(xy_normalized(:,1));
XYZ_hom = makehomogeneous(XYZ);
disp(size(XYZ_normalized));
disp(size(P));
XYZ_P = P*XYZ_hom;
disp(XYZ_P );
for i = 1:n
    error = error + norm(xy(:,i) - XYZ_P([1,2],i)/XYZ_P(3,i))^2;
end
error = sqrt(error/n);
disp('error DLT');
disp(error);
P = [P; 0,0,0,1];
end