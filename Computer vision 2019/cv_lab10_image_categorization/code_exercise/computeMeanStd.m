function [mu, sigma] = computeMeanStd(vBoW)
%      vBoW = [1,2,3,4,5;
%              2,3,3,3,1;
%              1,2,3,5,1];
    mu = mean(vBoW);
    sigma = std(vBoW);
end