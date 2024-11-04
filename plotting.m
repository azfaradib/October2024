[num,txt,raw] = xlsread('output_table.xlsx');
features=num(1:end,9:end);  

for i=1:size(features,1)
data =features(i,:); % Generate 1000 random samples from a normal distribution

% Create a histogram
[counts, edges] = histcounts(data, 'Normalization', 'pdf');

% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Calculate mean (μ)
mean(i)= sum(bin_centers .* counts) * (edges(2) - edges(1));

% Calculate standard deviation (σ)
variance = sum(((bin_centers - mean_value).^2) .* counts) * (edges(2) - edges(1));
std_dev(i) = sqrt(variance);
end

`