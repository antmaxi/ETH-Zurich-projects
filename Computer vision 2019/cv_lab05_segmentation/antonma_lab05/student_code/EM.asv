function [map, cluster] = EM(img)

% Number of Gaussians
K = 3;
% use function generate_mu to initialize mus
mu = generate_mu(img, K);
% use function generate_cov to initialize covariances
var = generate_cov(img, K);
% uniform mixture initially
alpha = 1/K*ones(1,K);


dim = size(img);
% aligning pixels in line
len = dim(1)*dim(2);
X = double(reshape(img, [len, 3]));

% iterate between maximization and expectation

% parameters
eps = 1e-4;
max_time_minutes = 20;
max_iter = 200;
iter = -0.5;
tic
log_like_old = -inf;
while true
    iter = iter + 1;
    P = expectation(mu,var,alpha,X);
    [mu, var, alpha] = maximization(P, X);
    time = toc/60;
    % create log-likelyhood value for stop criterion
    inverse = zeros(3,3,K);
    sqr_det = zeros(K,1);
    log_like = 0;
    for j = 1:K
        inverse(:,:,j) = inv(var(:,:,j));
        sqr_det(j) = det(var(:,:,j))^0.5;
    end
    for i = 1:len
        sum = 0;
        for j = 1:K
           pdf = alpha(j)*(2*pi)^(-3/2)/sqr_det(j)*exp(-0.5*(X(i,:)'-mu(:,j))'*inverse(:,:,j)*(X(i,:)'-mu(:,j)));
           sum = sum + pdf;
        end
        log_like = log_like + log(sum);
    end
    disp(log_like);
    if abs(log_like-log_like_old) < eps || iter > max_iter || time > max_time_minutes
        break;
    end
    log_like_old = log_like;
end
% output final parameters
disp('mu');
disp(mu);
disp('var');
disp(var);
disp('alpha');
disp(alpha);
disp('iter');
disp(iter);
disp('log_like');
disp(log_like_old);
map = zeros(dim(1), dim(2));
prob = zeros(dim(1), dim(2), K);
for i = 1:K
    prob(:,:,i) = reshape(P(:,i), dim(1), dim(2), 1);
end
% get peak with the maximal probabilty for each pixel
[M, map(:,:)] = max(prob(:,:,:), [], 3);

% assign colors for clusters (very contrast in order to see difference
% better)
if K == 3
    cluster = [1,0,0;
               0,1,0;
               0,0,1];
elseif K == 4
    cluster = [1,0,0;
               0,1,0;
               0,0,1;
               1,1,1];
elseif K == 5
    cluster = [1,0,0;
               0,1,0;
               0,0,1;
               1,1,1;
               0,0,0];
else
    cluster = rand(K,3);
end
end