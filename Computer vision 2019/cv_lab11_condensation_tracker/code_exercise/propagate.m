function particles = propagate(particles,sizeFrame,params)

% initialize model matrix
if (params.model == 0) %(0 = no motion, 1 = constant velocity)
    A = [1, 0;
         0, 1];
     
     coord_noise = normrnd(0, params.sigma_position, size(particles, 1), 2);
     particles = (A*particles')' + coord_noise;
     %disp(particles);
else
    A = [1, 0, 1, 0;
         0, 1, 0, 1;
         0, 0, 1, 0;
         0, 0, 0, 1];
     coord_noise = normrnd(0, params.sigma_position, size(particles, 1), 2);
     velo_noise =  normrnd(0, params.sigma_velocity, size(particles, 1), 2);  
     particles = (A*particles')' + [coord_noise velo_noise];
end

particles(:,1) = min(max(particles(:,1), 1), sizeFrame(2));
particles(:,2) = min(max(particles(:,2), 1), sizeFrame(1));
%disp(particles);

