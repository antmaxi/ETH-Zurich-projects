%% 
% Before starting, you should include the VLfeat (http://www.vlfeat.org/)
% and GCMex packages (https://github.com/shaibagon/GCMex) in your path:
% - for vlfeat under VLF_ROOT/vlfeat-0.9.21/toolbox/ you run vl_setup, 
% which does the job for you,
% - for GCMex under GCM_ROOT/ you run compile_gc, and then do 
% addpath('GCM_ROOT').
% Should you have any problems compiling them under Linux, Windows, Mac, 
% please refer to the corresponding websites for further instructions.

%%
%don't forget to initialize VLFeat
run(strcat(pwd, '\vlfeat-0.9.21\toolbox\vl_setup'));
%run(strcat(pwd, '\GCMex-master\GCMex-master'));
addpath(strcat(pwd, '\GCMex-master\GCMex-master'));
% Rectify images
imgNameL = 'images/0018.png';
imgNameR = 'images/0019.png';
camNameL = 'images/0018.camera';
camNameR = 'images/0019.camera';

scale = 0.5^2; % try this scale first
% scale = 0.5^3; % if it takes too long for GraphCut, switch to this one

imgL = imresize(imread(imgNameL), scale);
imgR = imresize(imread(imgNameR), scale);

figure(1);
subplot(121); imshow(imgL);
subplot(122); imshow(imgR);

[K R C] = readCamera(camNameL);
PL = K * [R, -R*C];
[K R C] = readCamera(camNameR);
PR = K * [R, -R*C];

[imgRectL, imgRectR, Hleft, Hright, maskL, maskR] = ...
    getRectifiedImages(imgL, imgR);

figure(2);
subplot(121); imshow(imgRectL);
subplot(122); imshow(imgRectR);

close all;
% Create an 15x15 square structuring element.
se = strel('square', 15);
maskL = imerode(maskL, se);
maskR = imerode(maskR, se);
%%
% Set disparity range
% (exercise 5.3)
% you may use the following two lines
%  [x1s, x2s] = getClickedPoints(imgRectL, imgRectR); 
%  close all;
% to get a good guess
% 
img1 = im2double(rgb2gray(imgL));
img2 = im2double(rgb2gray(imgR));
% Leave this exercise for the end, and for now try these fixed ranges
%extract SIFT features and match

% get disparities of points of interest from RANSAC on extracted by SIFT
% matches. Approximate range of disparities as maximal shift among this
% matching points along horizontal axis multiplied by 1.5 in order to make
% it more robust against existence of areas with no points of interest but
% bigger dispairity value
needed_number_of_inliers = 50;
number_of_inliers = 0;
threshold = scale*1e-6;

% iterate increasing threshold while we don't have enough matches
% while number_of_inliers < needed_number_of_inliers
%     [x1s, da] = extractSIFT(img1);
%     [x2s, db] = extractSIFT(img2);
%     da = single(da);
%     db = single(db);
%     [matches, scores] = vl_ubcmatch(da, db);
% 
%     x1s = [x1s(1:2,matches(1,:)); ones(1, size(matches,2))];
%     x2s = [x2s(1:2,matches(2,:)); ones(1, size(matches,2))];
% 
%     clf
%     
%     [F, inliers] = fundamentalMatrixRANSAC(x1s, x2s, threshold);
%     showFeatureMatches(img1, x1s(1:2, inliers), img2, x2s(1:2, inliers), 1);
%     max_disp = ceil(max(abs(x1s(2, inliers) - x2s(2, inliers)))*2);
%     dispRange = -max_disp:max_disp;
%     
%     number_of_inliers = size(inliers, 2);
%     disp("max dispairity");
%     disp(max_disp);
%     disp("number of inliers");
%     disp(number_of_inliers);
%     disp("___________________");
%     threshold = threshold * 10;
% end
dispRange = -40:40;
%%
% Compute disparities, winner-takes-all
% (exercise 5.1)
tic
hsize = 10;
%iterate over size of filter to get results with different parameters
for hsize = [3,7,12]
dispStereoL = ...
    stereoDisparity(rgb2gray(imgRectL), rgb2gray(imgRectR), dispRange, hsize);
dispStereoR = ...
    stereoDisparity(rgb2gray(imgRectR), rgb2gray(imgRectL), dispRange, hsize);
%disp(dispStereoR);
figure(1);
subplot(121); imshow(dispStereoL, [dispRange(1) dispRange(end)]);
subplot(122); imshow(dispStereoR, [dispRange(1) dispRange(end)]);
saveas(1, strcat('results\SD_', num2str(hsize), '.png'));
 thresh = 8;

maskLRcheck = leftRightCheck(dispStereoL, dispStereoR, thresh);
maskRLcheck = leftRightCheck(dispStereoR, dispStereoL, thresh);

maskStereoL = double(maskL).*maskLRcheck;
maskStereoR = double(maskR).*maskRLcheck;

figure(2);
subplot(121); imshow(maskStereoL);
subplot(122); imshow(maskStereoR);
saveas(1, strcat('results\SD_mask_', num2str(hsize), '.png'));
% close all;
toc
%%
% Compute disparities using graphcut
% (exercise 5.2)
 tic
% Labels = ...
%     gcDisparity(rgb2gray(imgRectL), rgb2gray(imgRectR), dispRange);
% dispsGCL = double(Labels + dispRange(1));
% Labels = ...
%     gcDisparity(rgb2gray(imgRectR), rgb2gray(imgRectL), dispRange);
% dispsGCR = double(Labels + dispRange(1));
% 
% figure(1);
% subplot(121); imshow(dispsGCL, [dispRange(1)+1+dispRange(end) dispRange(end)+1+dispRange(end)]);
% subplot(122); imshow(dispsGCR, [dispRange(1)+1+dispRange(end) dispRange(end)+1+dispRange(end)]);
% 
% maskLRcheck = leftRightCheck(dispsGCL, dispsGCR, thresh);
% maskRLcheck = leftRightCheck(dispsGCR, dispsGCL, thresh);
% 
% maskGCL = double(maskL).*maskLRcheck;
% maskGCR = double(maskR).*maskRLcheck;
% 
% figure(2);
% subplot(121); imshow(maskGCL);
% subplot(122); imshow(maskGCR);
% % % close all;
 toc
%%
% Using the following code, visualize your results from 5.1 and 5.2 and 
% include them in your report 
dispStereoL = double(dispStereoL);
dispStereoR = double(dispStereoR);
%  dispsGCL = double(dispsGCL);
%  dispsGCR = double(dispsGCR);
 
S = [scale 0 0; 0 scale 0; 0 0 1];
 
% For each pixel (x,y), compute the corresponding 3D point 
% use S for computing the rescaled points with the projection 
% % matrices PL PR
% [coords ~] = ...
%     generatePointCloudFromDisps(dispsGCL, Hleft, Hright, S*PL, S*PR);
%%same for other winner-takes-all
[coords2 ~] = ...
   generatePointCloudFromDisps(dispStereoL, Hleft, Hright, S*PL, S*PR);

imwrite(imgRectL, 'imgRectL.png');
imwrite(imgRectR, 'imgRectR.png');

% Use meshlab to open generated textured model, i.e. modelGC.obj
% generateObjFile('modelGC', 'imgRectL.png', ...
%     coords, maskGCL.*maskGCR);
% same for other winner-takes-all, i.e. modelStereo.obj
generateObjFile(strcat('new_modelStereo_', num2str(hsize)), 'imgRectL.png', ...
   coords2, maskStereoL.*maskStereoR);
 end