function particles_w = observe(particles,frame,H,W,hist_bin, ...
hist_target,sigma_observe)
xm = size(frame, 2);
ym = size(frame, 1);
particles_w = zeros(size(particles, 1));
for i = 1:size(particles, 1)
    xMin = min(max(1, round(particles(i,1) - W/2)), xm);
    xMax = min(max(1, round(particles(i,1) + W/2)), xm);
    yMin = min(max(1, round(particles(i,2) - H/2)), ym);
    yMax = min(max(1, round(particles(i,2) + H/2)), ym);

    hist = color_histogram(xMin,yMin,xMax,yMax,frame,hist_bin);

    chi = chi2_cost(hist, hist_target); 

    particles_w(i) = normpdf(chi,0,sigma_observe);
end
% normalize to have sum 1
particles_w = particles_w/sum(particles_w);
disp('sum_of_weights');
disp(sum(particles_w));