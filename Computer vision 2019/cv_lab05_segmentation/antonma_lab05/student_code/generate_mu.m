% Generate initial values for mu - uniformly random in min-max range
% K is the number of segments

function mu = generate_mu(imglab, K)
% get maximal and minimal values of pixels
maximal = double(max(max(imglab)));
minimal = double(min(min(imglab)));
mu = zeros(3, K);
for i = 1:3
    mu(i,:) = rand(1,K)*(maximal(1,1,i) - minimal(1,1,i)) + minimal(1,1,i);
end

end