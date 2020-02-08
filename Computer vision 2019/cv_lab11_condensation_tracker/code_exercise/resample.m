function [particles_new, particles_w_new] = resample(particles,particles_w)

% wheel sampling

n = size(particles, 1);
% Twice maximal weight
w_max2 = 2*max(particles_w);
% Starting uniformly random index of particle
index = randi([1 n],1,1);
% Array of uniformly random steps as a ratios to w_max2
rands = rand(n);
% Starting offset
beta = 0;
particles_new = particles;
particles_w_new = particles_w;
% Looping in "resampling wheel"
for i = 1:n
    beta = beta + rands(i)*w_max2;
    % Get closest index less than inside which we are at the moment
    while beta > particles_w(index)
        beta = beta - particles_w(index);
        index = mod(index, n) + 1;
    end
    particles_new(i) = particles(index);
    particles_w_new(i) = particles_w(index);
end
particles_w_new = particles_w_new/sum(particles_w_new);