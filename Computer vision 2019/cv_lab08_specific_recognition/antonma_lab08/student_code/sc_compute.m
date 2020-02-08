function [ShapeDescriptors1] = sc_compute(X,nbBins_theta,nbBins_r,smallest_r,biggest_r, tangent_flag, local_flag, deleted)

ShapeDescriptors1 = zeros(size(X,2), nbBins_theta, nbBins_r);
X = X(:, 1:size(X,2)-deleted);
mean_dist2 = mean(sqrt(dist2(X,X)), 'all')^2;
smallest_r2 = smallest_r * smallest_r;
biggest_r2 = biggest_r * biggest_r;
n = size(X,2);
disp(n);
disp(max(X(1,:)));
center = [(max(X(1,:)) + min(X(1,:)))/2, (max(X(2,:)) + min(X(2,:))/2)];
i = 0;
for p1 = X
    i = i + 1;
    min_dist_1 = Inf;
    min_dist_2 = Inf;
    ind_1 = 0;
    ind_2 = 0;
    j = 0;
    theta_array = zeros(1, n);
    r_array_ind = zeros(1, n) - 1;
    
    for p2 = X
        j = j + 1;
        % check that different point
        dp = p2 - p1;
        dist_sqr = (dp(1)^2+dp(2)^2)/mean_dist2;
        if (dist_sqr > 1e-8)   
            % get two closest points (may be used for angle calculation)
            if (dist_sqr <= min_dist_1)
                min_dist_2 = min_dist_1;
                ind_2 = ind_1;
                min_dist_1 = dist_sqr;
                ind_1 = j;         
            elseif ((dist_sqr > min_dist_1) && (dist_sqr < min_dist_2))
                min_dist_2 = dist_sqr;
                ind_2 = j;
            end
            %disp(dist_sqr);
            if ((dist_sqr >= smallest_r2) && (dist_sqr <= biggest_r2))
                r_array_ind(j) = 1 + floor(real(nbBins_r*(log(dist_sqr) - log(smallest_r2))/...
                                            (log(biggest_r2) - log(smallest_r2))));
                % check verticals
                if dp(1) == 0
                    if dp(2) > 0
                        theta = pi/2;
                    else
                        theta = -pi/2;
                    end
                else
                    if	local_flag
                        theta = atan(dp(2)/dp(1)); % if angles are measured relatively to local arrangement
                    else
                        
                        theta = atan((p2(2) - center(2))/(p2(1) - center(1))) -  atan((p1(2) - center(2))/(p1(1) - center(1))); % if reference point for angles is (0,0)
                    end
                end
                theta_array(j) =  real(wrapTo2Pi(real(theta)));%real(min(wrapTo2Pi(real(theta)), 2*pi - wrapTo2Pi(real(theta))));   
            end
            
        end
        
    end
    if (tangent_flag)
        theta_array = theta_array - theta_array(ind_1); % subtract tangent angle, 
    end                                                 % if angles are measured relatively to local arrangement

    theta_index = ones(1, n) + floor(theta_array/(2*pi+1e-8)*nbBins_theta);
    for j = 1:n
        if ((i~=j) && (r_array_ind(j) ~= -1))
            ShapeDescriptors1(i,  theta_index(j), r_array_ind(j)) = 1 + ShapeDescriptors1(i, theta_index(j), r_array_ind(j));
        end
    end
    
end
   
end
