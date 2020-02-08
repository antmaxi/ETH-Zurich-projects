% Compute the fundamental matrix using the eight point algorithm and RANSAC
% Input
%   x1, x2 	  Point correspondences 3xN
%   threshold     RANSAC threshold
%
% Output
%   best_inliers  Boolean inlier mask for the best model
%   best_F        Best fundamental matrix
%
function [best_inliers, best_F] = ransac8pF(x1, x2, threshold)

iter = 1000;

num_pts = size(x1, 2); % Total number of correspondences
best_num_inliers = 0; best_inliers = [];
best_F = 0;
M = iter;
for i=1:iter
    % Randomly select 8 points and estimate the fundamental matrix using these.
     [x1s,idx] = datasample(x1, 8, 2);
     x2s = x2(:, idx);
     [Fh, F] = fundamentalMatrix(x1s, x2s)
    % Compute the Sampson error.
    d1 = distPointsLines(x2, Fh*x1);
    d2 = distPointsLines(x1, Fh'*x2);   
    dist = d1 + d2;
    % Compute the inliers with errors smaller than the threshold.
    curr_inliers = abs(dist) < threshold;
    num_inliers = sum(curr_inliers);    
    % Update the number of inliers and fitting model if the current model
    % is better.
    if num_inliers > best_num_inliers
       best_num_inliers = num_inliers;
       best_F = Fh;
       best_inliers = curr_inliers;
       mean_samp = mean(dist(best_inliers));
    end
    r = best_num_inliers/num_pts;
    % Check if it's high probable that already was all inliers among 8
    % points
    % "adaptive" changes from adaptive RANSAC to usual one
    adaptive = true;
    if (1 - (1 - r^8)^i > 0.99) && adaptive
        M = i;
        break
    end    
end

disp("Total points");
disp(num_pts);
disp("Best number of inliers");
disp(best_num_inliers);
disp("Mean Sampson distance of inliers");
disp(mean_samp);
disp("Their ratio");
disp(best_num_inliers/mean_samp);
disp("Used threshold");
disp(threshold);
disp('M');
disp(M);
end


