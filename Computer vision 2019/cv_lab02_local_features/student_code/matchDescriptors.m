% Match descriptors.
%
% Input:
%   descr1        - k x q descriptor of first image
%   descr2        - k x q' descriptor of second image
%   matching      - matching type ('one-way', 'mutual', 'ratio')
%   
% Output:
%   matches       - 2 x m matrix storing the indices of the matching
%                   descriptors
function matches = matchDescriptors(descr1, descr2, matching)
    distances = ssd(descr1, descr2);
    %disp(distances);
    matches = [];
    if strcmp(matching, 'one-way')
        disp('one-way');
        [n, m] = size(distances);
        for i = 1:n
            [M,I] = min(distances(i,:));
            matches = [matches; i, I];
        end
        matches = matches.';
    elseif strcmp(matching, 'mutual')
        disp('mutual');
        [n, m] = size(distances);
        % Array of flags for descriptors from second image which are
        % already with pair from first (then 1 in array, otherwise 0)
        im = zeros(1, m);
        matches = [];
        % Get matches for the first image's corners
        for i = 1:n
            [M,I] = min(distances(i,:));
            if im(I) ~= 1
                [M1,I1] = min(distances(:,I));
                if I1 == i
                    im(I) = 1;
                    matches = [matches; i, I];
                end
           
            end
        end
        matches = matches.';
    elseif strcmp(matching, 'ratio')
        disp('ratio');
        % Ratio of minimal accepted difference in distances between 1st and
        % 2nd closest patches
        ratio = 0.5;
        matches = [];
        [n, m] = size(distances);
        disp('resulting ratios:');
        for i = 1:n
            [M,I] = mink(distances(i,:), 2);
            % Process case of existing only one patch on the second image
            size_I = size(I);
            if size_I(2) < 2
                matches = [matches; i, I];
            else
                % Check if patches are identical or 2nd is much closer 
                % than 1st patch => then add pair directly,
                % otherwise check ratio and add 1ss, if it's much closer than
                % 2nd
                r = distances(i,I(1))/distances(i,I(2));
                disp(r);
                if (distances(i,I(2)) - 0 < 10^(-8)) || (r > 1/ratio)
                    matches = [matches; i, I(2)];
                elseif (r < ratio)
                    matches = [matches; i, I(1)];
                end
            end
        end
        matches = matches.';
    else
        error('Unknown matching type.');
    end
end

function distances = ssd(descr1, descr2)
    % Numbers of descriptors in images
    s1 = size(descr1);
    n1 = s1(1);
    s2 = size(descr2);
    n2 = s2(1);
    % Prepare result matrix
    distances = zeros(n1, n2);
    % Iteration over all pairs of descriptors-patches
    for i = 1:n1
        patch1 = reshape(descr1(i,:,:), [1, 81]);
        for j = 1:n2
            patch2 = reshape(descr2(j,:,:), [1, 81]);
            % Distance as a sum of squared euclidian norms
            distances(i, j) = sum(pdist2(patch1, patch2, 'squaredeuclidean'), 'all');
        end
    end
end