function vBoW = create_bow_histograms(nameDir, vCenters)

  vImgNames = dir(fullfile(nameDir,'*.png'));
  nImgs = length(vImgNames);  
  vBoW  = zeros(0, size(vCenters,1));
  
  cellWidth = 4;
  cellHeight = 4;
  nPointsX = 10;
  nPointsY = 10;
  border = 8;
  
  % Extract features for all images in the given directory
  for i=1:nImgs 
    disp(strcat('  Processing image ', num2str(i),'...'));
    
    % load the image
    im1 = imread(fullfile(nameDir,vImgNames(i).name));
    % resize if needed
    scale = 1;
    im1 = imresize(im1, scale);
    img = double(rgb2gray(im1));
    % Collect local feature points for each image
    % and compute a descriptor for each local feature point
    vPoints = grid_points(img,nPointsX,nPointsY,border);
    [descriptors,~] = descriptors_hog(img,vPoints,cellWidth,cellHeight);
    
    % Create a BoW activation histogram for this image
    histo = bow_histogram(descriptors, vCenters);
    vBoW(i,:) = histo;
  end
    
end