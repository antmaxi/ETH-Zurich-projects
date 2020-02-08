function diffs = diffsGC(img1, img2, dispRange)

hsize = 3;
% get data costs for graph cut
diffs = zeros(size(img1,1), size(img1,2), size(dispRange,2));
for d = dispRange
    img1_shift = shiftImage(img1, d);
    ssd = (img1_shift(:,:)-img2(:,:)).^2;
    filter = fspecial('average', hsize);
    filtered = conv2(ssd, filter, 'same');
    diffs(:, :, d + max(dispRange) + 1) = filtered;
end
end