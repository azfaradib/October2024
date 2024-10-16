% Parameters for the two Gaussian curves
function [area_between] = findarea(a,b,c,d)
% a = 0;    % Mean of the first curve
% c = 1;    % Standard deviation of the first curve
% b = 2;    % Mean of the second curve
% d = 1.5;  % Standard deviation of the second curve

% Gaussian function definition
gaussian = @(x, mean, stdDev) (1 / (stdDev * sqrt(2 * pi))) * exp(-0.5 * ((x - mean) / stdDev) .^ 2);

% Function to find intersections
intersection_func = @(x) gaussian(x, a, c) - gaussian(x, b, d);

% Initial guesses for the intersection points
initial_guesses = linspace(-5, 5, 100);
intersection_points = [];

% Find intersection points
for guess = initial_guesses
    x = fzero(intersection_func, guess);
    if isempty(find(abs(intersection_points - x) < 1e-5, 1))  % Check for duplicates
        intersection_points = [intersection_points; x]; %#ok<AGROW>
    end
end

% Unique intersection points
intersection_points = unique(round(intersection_points, 5));

% Define the range for x values for integration
x_range = linspace(min(intersection_points), max(intersection_points), 1000);
y1 = gaussian(x_range, a, c);
y2 = gaussian(x_range, b, d);

% Calculate the area between the curves
area_between = trapz(x_range, abs(y1 - y2));
end


% % Plotting the curves
% figure;
% plot(x_range, y1, 'b', 'LineWidth', 2); hold on;
% plot(x_range, y2, 'r', 'LineWidth', 2);
% 
% % Fill the area between the curves
% fill([x_range, fliplr(x_range)], [y1, fliplr(y2)], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% title('Area Between Two Gaussian Curves');
% xlabel('x');
% ylabel('Probability Density');
% legend('Curve 1', 'Curve 2', 'Area Between Curves');
% grid on;
% 
% % Display the area
% disp('Area between the curves:');
% disp(area_between);