% Normalization of 2d-pts
% Inputs: 
%           xs = 2d points
% Outputs:
%           nxs = normalized points
%           T = 3x3 normalization matrix
%               (s.t. nx=T*x when x is in homogenous coords)
function [nxs, T] = normalizePoints2d(xs)
    % make mean = 0
    mu = mean(xs, 2);
    xs_shifted = xs - mu;
    % calculate axis rescaling factors
    T = zeros(3,3);
    for i = 1:2
        if norm(xs_shifted(i,:)) == 0
            scale = 0;
        else
            scale = sqrt(size(xs, 2))/norm(xs_shifted(i,:));
        end
        T(i, i) = scale;
        T(i, 3) = - mu(i)*scale;
    end
    T(3,3) = 1;
    % Get normalized points
    nxs = T*xs;
end
