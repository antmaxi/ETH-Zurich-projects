function [descriptors,patches] = descriptors_hog(img,vPoints,cellWidth,cellHeight)

    nBins = 8; % assume even
    nCellsW = 4; % number of cells, hard coded so that descriptor dimension is 128
    nCellsH = 4; 

    w = cellWidth; % set cell dimensions
    h = cellHeight;   

    pw = w*nCellsW; % patch dimensions
    ph = h*nCellsH; % patch dimensions

    descriptors = zeros(0,nBins*nCellsW*nCellsH); % one histogram for each of the 16 cells
    patches = zeros(0,pw*ph); % image patches stored in rows    
    
    [grad_x,grad_y] = gradient(img);    
    Gdir = (atan2(grad_y, grad_x)); 
    
    x_shifts = -nCellsH/2 : nCellsH/2 - 1;
    y_shifts = -nCellsW/2 : nCellsW/2 - 1;
    % shifts of left bottom corners of cells in nCellsW*nCellsH x 2 array
    [m, n] = ndgrid(x_shifts*nCellsH, y_shifts*nCellsW);
    cell_shifts = [m(:),n(:)];
    % shifts of pixels inside cell in cellHeight*cellWidth x 2 array
    [m, n] = ndgrid(0:cellHeight-1, 0:cellWidth-1);
    pixel_shifts = [m(:),n(:)];

    for i = 1:size(vPoints,1) % for all local feature points
        % coordinates of left bottom corners of cells
        cells = repmat(vPoints(i,:),nCellsW*nCellsH,1) + cell_shifts;
        %disp(cells);
        % store histograms for pixelsin cell
        cell_descriptor = zeros(0, nBins);
        % number of checked cells
        for cell_i = 1:size(cells, 1)
            pixels = repmat(cells(cell_i,:), cellHeight*cellWidth, 1)  + pixel_shifts;
            hist = zeros(1, nBins);
            %disp(pixels);
            for pixel_i = 1:size(pixels, 1)
                direction = floor(Gdir(pixels(pixel_i,1), pixels(pixel_i,2))/(2*pi/nBins + 1e-10)) + nBins/2 + 1;
                hist(direction) = 1 + hist(direction);
            end
            cell_descriptor(cell_i, :) = hist;
        end
        descriptors(i, :) = reshape(cell_descriptor, [1, nBins*nCellsW*nCellsH]);
        
        x_shifts = - nCellsH*cellHeight/2 : nCellsH*cellHeight/2 - 1;
        y_shifts = - nCellsW*cellWidth/2 : nCellsW*cellWidth/2 - 1;
        [m, n] = ndgrid(x_shifts, y_shifts);
        patch_shifts = [m(:),n(:)] + repmat(vPoints(i,:), pw*ph, 1);
        x = patch_shifts(:,1);
        y = patch_shifts(:,2);
        for j = 1:pw*ph
            patches(i, j) = img(x(j), y(j));
        end
    end % for all local feature points
end
