% Extract Harris corners.
%
% Input:
%   img           - n x m gray scale image
%   sigma         - smoothing Gaussian sigma
%                   suggested values: .5, 1, 2
%   k             - Harris response function constant
%                   suggested interval: [4e-2, 6e-2]
%   thresh        - scalar value to threshold corner strength
%                   suggested interval: [1e-6, 1e-4]
%   
% Output:
%   corners       - 2 x q matrix storing the keypoint positions
%   C             - n x m gray scale image storing the corner strength
function [corners, C] = extractHarris(img, sigma, k, thresh)
 grad_x = [-1/2, 0, 1/2];
 grad_y = [1/2, 0, -1/2].';
 I_x = conv2(grad_x, img);
 I_y = conv2(grad_y, img);
 %heatmap(I_x)
 %heatmap(I_y);
 % Gaussian kernel
 w = imgaussfilt([0,0,0;
                  0,1,0;
                  0,0,0], sigma);
 %disp(w);
 [n_x, n_y] = size(img);
 % Prepare matrix for Harris response function C(i,j)
 C = zeros(n_x, n_y);
 % Calculate C(i,j)
 for i = 2:(n_x - 1)
     for j = 2:(n_y - 1)
         M = [0,0;
              0,0];
         for di = 1:3
              for dj = 1:3    
                I_x_loc = I_x(i + di - 2, j + dj - 2); 
                I_y_loc = I_y(i + di - 2, j + dj - 2);
                M = M + w(di,dj)*[I_x_loc*I_x_loc, I_x_loc*I_y_loc;
                                  I_x_loc*I_y_loc, I_y_loc*I_y_loc];
              end
         end
         e = eig(M);
         % C(i, j) = det(M)-k*Tr(M)^2
         C(i, j) = e(1)*e(2) - k*(e(1)+e(2))^2;
     end
 end

 % Find maximum in 3x3 neighborhoods of C
 corners = [];
 Max_reg = imregionalmax(C, 8);
 for i = 2:(n_x - 1)
     for j = 2:(n_y - 1)
         if (abs(C(i, j)) > thresh) && (Max_reg(i, j))
             corners = [corners; i, j];
         end
     end
 end
 corners = corners.';
end