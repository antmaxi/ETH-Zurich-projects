function meanState = estimate(particles,particles_w)
% particles = [1,1;  
%              0,0;
%              2,3];
% particles_w = [0.2, 0.2, 0.6];

% disp(particles);
% disp('particles_w');
% disp(particles_w);
% disp('sizeParticles');
% disp(size(particles));
% disp('particles');
% disp(particles(:,1));
% disp('sumparticles');
% disp(sum(particles(:,1).*particles_w(:)));
meanState = zeros(1, size(particles, 2));
for i = 1:size(particles, 2)
    meanState(i) = sum(particles(:,i).*particles_w(:));
end
%disp(particles);
% disp('meanState');
% %disp(size(meanState));
% disp(meanState);