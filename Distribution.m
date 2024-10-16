% Create a normally distributed (mu: 5, sigma: 3) random data set
x = normrnd(5, 3, 1e4, 1);

% Compute and plot results. The results are sorted by "Bayesian information
% criterion".
[D, PD] = allfitdist(x, 'PDF');