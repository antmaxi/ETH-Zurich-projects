function [T] = transform_matrix(xy)
% get mean, transformation matrix for
% non-homogenizated data
s = size(xy);

n = s(1, 2);
m = s(1, 1);
% 1. compute centroid
mu = mean(xy.');

% 2. shift the input points so that the centroid is at the origin
mu_rep = repmat(mu.', 1, n);
xy = xy - mu_rep;

% 3. compute scale

T = zeros(m, m);
for j = 1:m
    scale = 0;
    for i = 1:n
        scale = scale + xy(j, i)^2;
    end
    T(j,j) = sqrt(n/scale);     
end
% 4. create T and U transformation matrices (similarity transformation)
T = [T, -mu.'; zeros(1, m), 1];
% unkomment to check what would be without notmalization
%T = eye(m+1);

end