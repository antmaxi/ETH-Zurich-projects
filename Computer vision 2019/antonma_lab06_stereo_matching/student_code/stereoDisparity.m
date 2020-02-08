function best_d = stereoDisparity(img1, img2, dispRange, hsize)
% dispRange: range of possible disparity values
% --> not all values need to be checked


% best dispairities for pixels
best_d = zeros(size(img1));
% smallest ssd for pixels
best_ssd = Inf(size(img1));
for d = dispRange
    img1_shift = shiftImage(img1, d);
    ssd = (img1_shift(:,:)-img2(:,:)).^2;
    filter = fspecial('average', hsize);
    filtered = conv2(ssd, filter, 'same');
    mask = filtered < best_ssd; 
    best_d(mask) = d +1+dispRange(end);
    best_ssd(mask) = ssd(mask);
end
disp(ssd)
end