function [best_k, best_b] = ransacLine(data, iter, threshold)
% data: a 2xn dataset with n data points
% iter: the number of iterations
% threshold: the threshold of the distances between points and the fitting line

num_pts = size(data, 2); % Total number of points
best_num_inliers = 0;   % Best fitting line with largest number of inliers
best_k = 0; best_b = 0;     % parameters for best fitting line

for i=1:iter
    % Randomly select 2 points and fit line to these
    % Tip: Matlab command randperm / randsample is useful here
    k = randsample(num_pts, 2);
    xy1 = data(:, k(1));
    xy2 = data(:, k(2));
    % Model is y = k*x + b. You can ignore vertical lines for the purpose
    % of simplicity.
    % sample while line is vertical
    while xy1(1) - xy2(1) == 0
        sample = randsample(num_pts, 2);
        xy1 = data(:, sample(1));
        xy2 = data(:, sample(2));
    end
    k = (xy1(2) - xy2(2))/(xy1(1) - xy2(1));
    b = xy1(2) - k * xy1(1);
    % Compute the distances between all points with the fitting line
    dist = (k*data(1, :) + b - data(2, :)) / sqrt(1 + k * k);    
    % Compute the inliers with distances smaller than the threshold
    inliers = sum(abs(dist) < threshold);
    % Update the number of inliers and fitting model if the current model
    % is better.
    if inliers > best_num_inliers
       best_num_inliers = inliers;
       best_k = k; 
       best_b = b;
    end
end


end
