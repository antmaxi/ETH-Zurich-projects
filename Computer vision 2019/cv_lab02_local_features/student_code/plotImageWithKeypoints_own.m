% Plot image with keypoints saving file with parameters in the name.
%
% Input:
%   img        - n x m color image 
%   keypoints  - 2 x q matrix, holding keypoint positions
%   fig        - figure id
function plotImageWithKeypoints(img, keypoints, fig, sigma, k, thresh, type)
    figure(fig);
    imshow(img, []);
    hold on;
    plot(keypoints(2, :), keypoints(1, :), '+r');
    hold off;
    sig = num2str(sigma);
    thr = num2str(thresh);
    k1 = num2str(k);
    saveas(fig,strcat('C:\Users\Anton\Desktop\ETH_books\CV\Lab_Assignment_02-_Files\student_code\res2\', type, '-', ...
        sig(end), '-', thr(1), '-', k1(end)),'png');
end