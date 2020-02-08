function vPoints = grid_points(img,nPointsX,nPointsY,border)
%     nPointsX = 3;
%     nPointsY = 2;
%     border = 8;
%     x_dim = 101;
%     y_dim = 50;
    
    [x_dim, y_dim] = size(img);
    vPoints = zeros(nPointsX * nPointsY, 2);
    
    x = 1:nPointsX;
    y = 1:nPointsY;
    [X,Y] = meshgrid(x,y);
    
    x_step = (x_dim - border * 2 - 1) / (nPointsY - 1);
    y_step = (y_dim - border * 2 - 1) / (nPointsX - 1);
    
    % round in order to get pixels values, reshape to nPointsX*nPointsY x 2
    % images in Matlab are with inverted x-y coordrinate system (x|, y-)
    vPoints(:,1) = reshape(round(x_step*(X-1) + border + 1), [nPointsX * nPointsY, 1]);
    vPoints(:,2) = reshape(round(y_step*(Y-1) + border + 1), [nPointsX * nPointsY, 1]);
    
end
