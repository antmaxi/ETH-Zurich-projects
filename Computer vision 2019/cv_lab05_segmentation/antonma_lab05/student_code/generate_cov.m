% Generate initial values for the K
% covariance matrices

function cov = generate_cov(imglab, K)
cov = zeros(3, 3, K);
% get maximal and minimal values of pixels (i.e. scale of the space color
% in which we are working)
maximal = double(max(max(imglab)));
minimal = double(min(min(imglab)));
for i=1:K
    cov(:,:,i) =         [maximal(1,1,1)-minimal(1,1,1), 0, 0;
                          0, maximal(1,1,2)-minimal(1,1,2), 0;
                          0, 0, maximal(1,1,3)-minimal(1,1,3)];
end
end