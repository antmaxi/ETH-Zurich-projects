function label = bow_recognition_bayes(histogram, vBoWPos, vBoWNeg)


[muPos, sigmaPos] = computeMeanStd(vBoWPos);
[muNeg, sigmaNeg] = computeMeanStd(vBoWNeg);

% Calculating the probability of appearance each word in observed histogram
% according to normal distribution in each of the positive and negative bag
% of words
log_prob_pos = 0;
log_prob_neg = 0;
for i = 1:size(histogram,2)
    % check that not degenerate normal distribution
    if (muPos(i) ~= 0) || (sigmaPos(i) ~= 0)
        log_prob_pos = log_prob_pos + log(normpdf(histogram(i), muPos(i), sigmaPos(i)));
    end
    if (muNeg(i) ~= 0) || (sigmaNeg(i) ~= 0)
        log_prob_neg = log_prob_neg + log(normpdf(histogram(i), muNeg(i), sigmaNeg(i)));
    end
end

p_car = 0.5; % prior
prob_pos = p_car * exp(log_prob_pos);
prob_neg =(1 - p_car) * exp(log_prob_neg);

if (prob_pos > prob_neg)
    label = 1;
else
    label = 0;
end