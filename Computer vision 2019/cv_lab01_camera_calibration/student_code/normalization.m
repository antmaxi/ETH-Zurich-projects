%% !!! DO NOT CHANGE THE FUNCTION INTERFACE, OTHERWISE, YOU MAY GET 0 POINT !!! %%
% xy: size 2xn
% XYZ: size 3xn
% xy_normalized: 3xn
% XYZ_normalized: 4x
function [xy_normalized, XYZ_normalized, T, U] = normalization(xy, XYZ)
%data normalization
% 1. compute centroid
% 2. shift the input points so that the centroid is at the origin
% 3. compute scale
% 4. create T and U transformation matrices (similarity transformation)
[T] = transform_matrix(xy);

[U] = transform_matrix(XYZ);

% 5. normalize the points according to the transformations
xy = homogenization(xy);
XYZ = homogenization(XYZ);
xy_normalized = T*xy;
XYZ_normalized = U*XYZ;

end