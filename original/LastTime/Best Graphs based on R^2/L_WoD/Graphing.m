% Clear workspace, command window, and close all figures
clear; clc; close all;

% --- 1. Read and Filter Data (Same as before) ---
data = readtable('L_WoD(40).xlsx');
time = data{:,1};
force = data{:,2};
window_size = 51; % Smoothing window size
filtered_force = smoothdata(force, 'sgolay', window_size);

% --- 2. Split the Data (Same as before) ---
split_time = time(end) / 2;
[~, split_index] = min(abs(time - split_time));

time_part1 = time(1:split_index);
force_part1 = filtered_force(1:split_index);

time_part2 = time(split_index:end);
force_part2 = filtered_force(split_index:end);

% --- 3. Fit a Separate Function to Each Piece (Same as before) ---
degree_part1 = 5;
degree_part2 = 17;

p1 = polyfit(time_part1, force_part1, degree_part1);
fitted_force_part1 = polyval(p1, time_part1);

p2 = polyfit(time_part2, force_part2, degree_part2);
fitted_force_part2 = polyval(p2, time_part2);


% --- 4. Plot the Piecewise Fit (Same as before) ---
figure('Name', 'Piecewise Function Fit');

plot(time, filtered_force, 'bo', 'MarkerSize', 8, 'DisplayName', 'Experimental Data');
hold on;
plot(time_part1, fitted_force_part1, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
plot(time_part2, fitted_force_part2, 'r-', 'LineWidth', 2, 'HandleVisibility','off');
xline(split_time, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');


% --- 5. Customize the Plot ---
% *** This section has been updated ***
ylim([-2e-9, 10e-9]);
xlabel('Time (s)');
ylabel('Force (N)');
title('Piecewise Fit of Experimental Data');

% Set legend location to top-right and turn off grid lines
legend('show', 'Location', 'northeast'); %'northeast' corresponds to top-right
grid off; % This command removes the grid lines

hold off;

% --- 6. Save the Combined Results (Same as before) ---
combined_fit = [fitted_force_part1; fitted_force_part2(2:end)];
final_data = table(time, filtered_force, combined_fit, 'VariableNames', {'Time_s', 'Filtered_Force', 'Piecewise_Fitted_Force'});
writetable(final_data, 'piecewise_fitted_data.xlsx');

disp('New Excel file "piecewise_fitted_data.xlsx" has been created.');