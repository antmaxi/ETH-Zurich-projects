% Extract descriptors.
%
% Input:
%   img           - the gray scale image
%   keypoints     - detected keypoints in a 2 x q matrix
%   
% Output:
%   keypoints     - 2 x q' matrix
%   descriptors   - w x q' matrix, stores for each keypoint a
%                   descriptor. w is the size of the image patch,
%                   represented as vector
function [keypoints, descriptors] = extractDescriptors(img, keypoints)
    k = size(keypoints);
    k = k(1, 2);
    [n_x, n_y] = size(img);
    descriptors = zeros(1,9,9);
    number = 0;
    % Array to store indices of close to the image boundary keypoints
    bad_keypoints = [];
    for i = 1:k
        x = keypoints(1,i);
        y = keypoints(2,i);
        % Add patches
        if ((x > 4) && (y > 4) && ((n_x - x) > 3) && ((n_y - y) > 3))   
           	number = number + 1;
            descriptors(number,:,:) = img(x-4:x+4, y-4:y+4);
        else
            bad_keypoints = [bad_keypoints, i];
        end
    end
    % Deleting close to the boundary points
    for i = size(bad_keypoints):1
        keypoints(:,i) = [];
    end
end