function hist = color_histogram(xMin,yMin,xMax,yMax,frame,hist_bin)
%disp(xMin);

im = imcrop(frame, [xMin, yMin, xMax-xMin, yMax-yMin]);
hist = imhist(im, hist_bin^3);
hist = hist/sum(hist);

%     x_range = xMin : xMax;
%     y_range = yMin : yMax;
%     [m, n] = ndgrid(x_range, y_range);
%     pixels_shifts = [m(:),n(:)];
%     pixels = zeros(size(pixels_shifts, 1), 3);
% %     disp('aaa');
% %     disp(xMin);
% %     disp(xMax);
% %     disp(yMin);
% %     disp(yMax);
% %     disp('aaa');
%     %pixel_curr = frame([pixels_shifts(,2), pixels_shifts(i,1), :);
%     for i = 1:size(pixels_shifts, 1) 
%         pixel_curr = frame(pixels_shifts(i,2), pixels_shifts(i,1), :);
%         pixels(i,:) = pixel_curr(1,1,:);
%     end
    
    %disp(size(pixels));
    %disp(pixels);

%     hist = zeros(hist_bin^3);
%     ind1 = floor(pixels(:,1)./(256/hist_bin))+1;
%     ind2 = floor(pixels(:,2)./(256/hist_bin))+1;
%     ind3 = floor(pixels(:,3)./(256/hist_bin))+1;
%     hist(ind1 + ind2*hist_bin + ind3*hist_bin^2) = 1 ...
%                     + hist(ind1 + ind2*hist_bin + ind3*hist_bin^2);
% %     for i = 1:size(pixels_shifts, 1)
%         ind1 = floor(pixels(i,1)/(256/hist_bin))+1;
%         ind2 = floor(pixels(i,2)/(256/hist_bin))+1;
%         ind3 = floor(pixels(i,3)/(256/hist_bin))+1;
%         hist(ind1 + ind2*hist_bin + ind3*hist_bin^2) = 1 ...
%                     + hist(ind1 + ind2*hist_bin + ind3*hist_bin^2);
%     end
    