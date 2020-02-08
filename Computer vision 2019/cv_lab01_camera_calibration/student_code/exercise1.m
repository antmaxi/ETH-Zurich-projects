% Exervice 1
%
clear all;
close all;

IMG_NAME = 'image001.jpg';

%This function displays the calibration image and allows the user to click
%in the image to get the input points. Left click on the chessboard corners
%and type the 3D coordinates of the clicked points in to the input box that
%appears after the click. You can also zoom in to the image to get more
%precise coordinates. To finish use the right mouse button for the last
%point.
%You don't have to do this all the time, just store the resulting xy and
%XYZ matrices and use them as input for your algorithms.
[xy XYZ] = getpoints(IMG_NAME);

%Data that was used for report
%first set
%XYZ =  [0,0,0,7,4,7,0,0,0,3,2,6;
%        6,3,6,0,0,0,4,5,2,0,0,0;
%        0,4,9,9,3,0,1,4,7,7,3,3];
%xy = [1414.9, 1043.8,    1233.1,    232.8,    493.3,    68.1, 1159.3,    1215.9,    957.8,    608.8,    692.4,    254.9;
%    1174.4,    665.6,    19.2,    129.8,    788.5,    1167.0, 1014.6,    606.6,    412.5,    419.9,    808.2,    773.8];
%second set
%xy_normalized =
%   1.0e+03 *
%   1.0709    1.1839    0.6997    0.5105    0.6874    0.0681
%    1.0195    0.3805    0.2822    0.4838    0.8106    1.1744
%    0.0010    0.0010    0.0010    0.0010    0.0010    0.0010
%XYZ_normalized =
%     0     0     2     4     2     7
%     3     5     0     0     0     0
%     1     6     9     6     3     0
%     1     1     1     1     1     1
disp('Data');
disp('xy');
disp(xy);
disp('XYZ');
disp(XYZ);
% === Task 1 Data normalization ===
[xy_normalized, XYZ_normalized, T, U] = normalization(xy, XYZ);

% === Task 2 DLT algorithm ===

[P, K, R, t, error] = runDLT(xy, XYZ)
visualization_reprojected_points(xy, XYZ, P, IMG_NAME);
%building array of coordinates of corners
XYZ_all = [];
for j = 0:6
    for k = 0:9
        vec = [0,j,k].';
        XYZ_all = [XYZ_all vec];
    end
end
for i = 0:7
    for k = 0:9
        vec = [i,0,k].';
        XYZ_all = [XYZ_all vec];
    end
end
%visualizing these corners
visualization_3d_points (XYZ_all, P, IMG_NAME)
% === Task 3 Gold Standard algorithm ===
[P, K, R, t, error] = runGoldStandard(xy, XYZ);
visualization_reprojected_points(xy, XYZ, P, IMG_NAME);
visualization_3d_points (XYZ_all, P, IMG_NAME)