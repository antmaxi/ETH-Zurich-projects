function [map, peak] = meanshiftSeg(img)

dim = size(img);
len = dim(1)*dim(2);
% parameter - radius
r = 5;

% aligning pixels in line
img_lin = reshape(img, [len, 3]);
map = zeros(len,1);
map(1) = 1;
% coordinates of peaks
peak =  zeros(len,3);
% how many pixels are merged already in corresponding peak
votes = zeros(len,3);
peak(1,:) = findPeak(img_lin, 1, r);
% Resulting number of peaks
K = 1;
for i = 2:len
    new_peak = findPeak(img_lin, i, r);
    merged = false;
    for j = 1:i-1
        if norm(new_peak - peak(map(j),:)) < r/2
            merged = true;
            peak(map(j),:) = (new_peak*votes(map(j)) + peak(map(j),:))/(votes(map(j))+1);
            votes(map(j)) = votes(map(j))+1;
            map(i) = map(j);
            break;
        end
    end
    if ~merged
        K = K + 1;
        map(i) = K;
        peak(K, :) = new_peak;
        votes(K) = 1;
    end
end
% final check if there are new close peaks
% close = true;
% while close
%     close = false;
%     for i=1:len
%         for j = 1:i-1
%             if norm(peak(i) - peak(j) < r/2
%                 close = true;
%                 peak(map(j),:) = (peak(map(j),:)*votes(map(i)) + peak(map(j),:)*votes(map(j)))/...
%                     (votes(map(j))+(votes(map(i))));
%                 map(i) = map(j);
%                 votes(map(j)) = votes(map(j))+1;
%                 break;
%             end
%         end
%         if close == true
%             break;
%         end
%     end
% end

                
                
peak = double(peak(1:K,:));
map = reshape(map, [dim(1), dim(2)]);
end