% show feature matches between two images
%
% Input:
%   img1        - n x m color image 
%   corner1     - 2 x k matrix, holding keypoint coordinates of first image
%   img2        - n x m color image 
%   corner1     - 2 x k matrix, holding keypoint coordinates of second image
%   fig         - figure id
function showFeatureMatches(img1, corner1, img2, corner2, fig, out1, out2)
    [sx, sy, sz] = size(img1);
    img = [img1, img2];
    
    corner2 = corner2 + repmat([sy, 0]', [1, size(corner2, 2)]);
    
    figure(fig), imshow(img, []);    
    hold on, plot(corner1(1,:), corner1(2,:), '+r');
    hold on, plot(corner2(1,:), corner2(2,:), '+r');    
    hold on, plot([corner1(1,:); corner2(1,:)], [corner1(2,:); corner2(2,:)], 'b');    
    
   	out2 = out2 + repmat([sy, 0]', [1, size(out2, 2)]);
    
    hold on, plot(out1(1,:), out1(2,:), '+b');
    hold on, plot(out2(1,:), out2(2,:), '+b');    
    hold on, plot([out1(1,:); out2(1,:)], [out1(2,:); out2(2,:)], 'r');
    %saveas(fig,strcat(num2str(fig), '.png'));
    close
end