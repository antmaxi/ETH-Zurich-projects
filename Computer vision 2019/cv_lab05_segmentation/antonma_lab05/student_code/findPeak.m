function peak = findPeak(X, pixel, r)
total_pixels = size(X,1);
% threshold when moving of sphere stops
threshold = 0.01;
% Shift of the center of sphere by one step
shift = inf;
% Center of the sphere
center = double(X(pixel,:));
while shift > threshold
    % Sum of coord of points inside sphere
    sum = double(zeros(1,3));
    % Number of points inside the sphere
    num = 0;
    % Checking pixels if they are in r-radius area
    for index = 1:total_pixels
       x = double(X(index, :));
       z = x - center;
       dist = norm(z);   
       if dist < r
          num = num + 1;
          sum = sum + x;
       end
    end
    new_center = sum/num;
    shift = norm(center - new_center);
    center = new_center;
end    
peak = double(center);
end