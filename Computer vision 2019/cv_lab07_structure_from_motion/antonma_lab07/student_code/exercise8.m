% =========================================================================
% Exercise 8
% =========================================================================
addpath('dispairity');
% Initialize VLFeat (http://www.vlfeat.org/)
run(strcat(pwd, '\vlfeat-0.9.21\toolbox\vl_setup'));
%K Matrix for house images (approx.)
K = [  670.0000     0     393.000
         0       670.0000 275.000
         0          0        1];

%Load images
imgName1 = '../data/house.000.pgm';
imgName2 = '../data/house.004.pgm';

img1 = single(imread(imgName1));
img2 = single(imread(imgName2));

%extract SIFT features and match
[fa, da] = vl_sift(img1);
[fb, db] = vl_sift(img2);
%don't take features at the top of the image - only background
filter = fa(2,:) > 100;
fa = fa(:,find(filter));
da = da(:,find(filter));

[matches, scores] = vl_ubcmatch(da, db);
%showFeatureMatches(img1, fa(1:2, matches(1,:)), img2, fb(1:2, matches(2,:)), 20);
x1s = fa(1:2, matches(1,:));
x2s = fb(1:2, matches(2,:));
%% Compute essential matrix and projection matrices and triangulate matched points

%use 8-point ransac or 5-point ransac - compute (you can also optimize it to get best possible results)
%and decompose the essential matrix and create the projection matrices
t = 1e-3;
% x1s - 2xN
[F, inliers] = ransacfitfundmatrix(x1s, x2s, t);
x1_selected = x1s(1:2,inliers);
x2_selected = x2s(1:2,inliers);
x1_out = x1_selected(1:2, setdiff(1:size(x1_selected,2), inliers));
x2_out = x2_selected(1:2, setdiff(1:size(x2_selected,2), inliers));
showFeatureMatches(img1, x1_selected, img2, x2_selected, 20, x1_out, x2_out);

x1_sel_hom = makehomogeneous(x1_selected);
x2_sel_hom = makehomogeneous(x2_selected);
% calibrated with camera matrices
x1_calibrated = K\x1_sel_hom;
x2_calibrated = K\x2_sel_hom;
x1_calibrated = hnormalise(x1_calibrated);
x2_calibrated = hnormalise(x2_calibrated);
% get essential matrix using info about camera matrices
E = K'*F*K;
n = 4;
% cell of projection matrices
Ps = cell(n+1);
% cell of 3D-point in common coordinates
XS = cell(n+1);
Ps{n+1} = eye(4);
Ps{n} = decomposeE(E, x1_calibrated, x2_calibrated);
%triangulate the inlier matches with the computed projection matrix
[XS{n+1}, err] = linearTriangulation(Ps{n+1}, x1_calibrated, Ps{n}, x2_calibrated);


%%  SHOW EPIPOLES
x1s = makehomogeneous(x1_selected);
x2s = makehomogeneous(x2_selected);
% Get right and left null vectors of fundamental matrix and normalize them
% to obtain positions of epipoles
[U,S,V] = svd(F);
epipole = V(:, end);
epipole(1) = epipole(1)/epipole(3);
epipole(2) = epipole(2)/epipole(3);

% show clicked points
figure(1),clf, imshow(img1, []); hold on, plot(x1s(1,:), x1s(2,:), '*r'),
%plot epipole in img 1
plot(epipole(1), epipole(2), 'go');

epipole = U(:, end);
epipole(1) = epipole(1)/epipole(3);
epipole(2) = epipole(2)/epipole(3);
figure(2),clf, imshow(img2, []); hold on, plot(x2s(1,:), x2s(2,:), '*r'),
%plot epipole in img 2
plot(epipole(1), epipole(2), 'go');

% draw epipolar lines in img 1
figure(1)
for k = 1:size(x1s,2)
    drawEpipolarLines(F'*x2s(:,k), img1);
end

% draw epipolar lines in img 2
figure(2)
for k = 1:size(x2s,2)
    drawEpipolarLines(F*x1s(:,k), img2);
end

%% Add an addtional views of the scene 
t = 1e-2;
% function to add projection matrix and 3D-points from new view image
 for i = 1:3 
     imgName = strcat('../data/house.00', num2str(i), '.pgm');
     [Ps{i}, XS{i}, err] = addView(i, imgName, da, x1_selected, x1_calibrated, XS{n+1}, K, t, Ps{n+1}, matches, inliers, img1);
 end
%% Plot stuff


%use plot3 to plot the triangulated 3D points.
% the first view points in red
fig = 10;
figure(fig);
plot3(XS{n+1}(1,:,:), XS{n+1}(2,:,:), XS{n+1}(3,:,:), '.r',  XS{1}(1,:,:), XS{1}(2,:,:), XS{1}(3,:,:), '.b', ...
      XS{2}(1,:,:), XS{2}(2,:,:), XS{2}(3,:,:), '.k', XS{3}(1,:,:), XS{3}(2,:,:), XS{3}(3,:,:), '.g');
%draw cameras
drawCameras(Ps, fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DENSE RECONSTRUCTION
    
    imgName1 = '../data/house.002.pgm';
    imgName2 = '../data/house.003.pgm';
    
    scale = 0.5;
    [dispair, imgRectL, coords] = exercise6(imgName1, imgName2, Ps{2}(1:3,:), Ps{3}(1:3,:), K, scale);

%% function to add views
function [Ps, XS, err] = addView(i, imgName, da, x1_selected, x1_calibrated, XS, K, t, P0, matches0, inliers, img1)
    img3 = single(imread(imgName));
    [fc, dc] = vl_sift(img3);
    %match against the features from image 1 that where triangulated
    [matches, ~] = vl_ubcmatch(da(:,matches0(1,inliers)), dc);
    x1_selected_new = x1_selected(1:2, matches(1,:));
    x3s = fc(1:2, matches(2,:));
    XSs = XS(1:3, matches(1,:));

    %run 6-point ransac with DLT inside (projmatrix.m)
    % calibrated with camera matrices
    x3s = makehomogeneous(x3s);
    x3_calibrated = K\x3s;
    x3_calibrated = hnormalise(x3_calibrated);
    [Ps, inliers1] = ransacfitprojmatrix(x3_calibrated, XSs, t);
    %%%%%%%%%% USE FOLLOWING LINES OF CODE TO FIX ORIENTATION OF CAMERAS %%%%%%%%%%
    if (det(Ps) < 0 )
        Ps(1:3,1:3) = -Ps(1:3,1:3);
        Ps(1:3, 4) = -Ps(1:3, 4);
    end

    % x3_show 2xN
    x3_show = x3s(1:2, inliers1);
    x3_out = x3s(1:2, setdiff(1:size(x3s,2), inliers1));
    x1_out = x1_selected_new(1:2, setdiff(1:size(x3s,2), inliers1));
    disp(size(x1_out));
    showFeatureMatches(img1, x1_selected_new(1:2, inliers1), img3, x3_show, 20+i, x1_out, x3_out);

    % XS_selected 3xN
    x3_inliers = x3_calibrated(1:2, inliers1);

    x1_calibrated_selected = x1_calibrated(1:3, inliers1);
    x1_calibrated_selected = hnormalise(x1_calibrated_selected);
    x3_sel_hom = makehomogeneous(x3_inliers);


    %triangulate the inlier matches with the computed projection matrix
    [XS, err] = linearTriangulation(P0, x1_calibrated_selected, Ps, x3_sel_hom);
end
