A = [5, 3, 8; 2, 7, 4; 6, 1, 9];

% Define the threshold value
threshold = 4;

% Find the indices of elements less than the threshold
[row, col] = find(A < threshold);

% Display the results
indices = [row, col];
disp('Indices of elements less than the threshold:');
disp(indices);