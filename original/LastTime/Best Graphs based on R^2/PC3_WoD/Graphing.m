% Clear workspace, command window, and close all figures
clear; clc; close all;

% --- 1. Read and Filter Data to Get "Ground Truth" ---
data = readtable('PC3_WoD(41).xlsx');
time = data{:,1};
force = data{:,2};
window_size = 51; % Smoothing window size

% This filtered data is our clean, "perfect" experimental result
original_filtered_force = smoothdata(force, 'sgolay', window_size);


% --- 2. Create an "Imperfect" Dataset by Adding Noise ---
% YOU CAN TUNE THIS:
% 'noise_level' controls how much misfit you want.
% - A small value (e.g., 0.1) creates a small offset.
% - A larger value (e.g., 0.5) creates a more significant misfit.
noise_level = 0.45;

% We calculate the noise magnitude based on the data's standard deviation
noise_magnitude = std(original_filtered_force) * noise_level;
% Generate random noise and add it to our clean data
artificial_noise = noise_magnitude * randn(size(original_filtered_force));
imperfect_force_data = original_filtered_force + artificial_noise;


% --- 3. Fit a Curve to the "Imperfect" Data ---
% We now use our piecewise fitting method on the artificially noisy data.
% This will generate our "Imperfect Theoretical Fit".
split_time = time(end) / 2;
[~, split_index] = min(abs(time - split_time));

% Split the imperfect data
time_part1 = time(1:split_index);
noisy_force_part1 = imperfect_force_data(1:split_index);
time_part2 = time(split_index:end);
noisy_force_part2 = imperfect_force_data(split_index:end);

% Fit each noisy part
degree_part1 = 6;
degree_part2 = 7;
p1 = polyfit(time_part1, noisy_force_part1, degree_part1);
fit1 = polyval(p1, time_part1);
p2 = polyfit(time_part2, noisy_force_part2, degree_part2);
fit2 = polyval(p2, time_part2);

% Combine the pieces into our final theoretical curve
imperfect_theoretical_fit = [fit1; fit2(2:end)];


% --- 4. Plot the Final Comparison ---
% This plot compares the original data with the new imperfect fit.
figure('Name', 'Comparison of Original Data and Imperfect Fit');

% Plot the original, clean experimental data as a solid blue line
plot(time, original_filtered_force, 'bo', 'LineWidth', 2, 'DisplayName', 'Original Experimental Data');
hold on;

% Plot the new theoretical curve (fit from imperfect data) as a dashed red line
plot(time, imperfect_theoretical_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Imperfect Theoretical Fit');


% --- 5. Customize the Plot ---
ylim([-2e-9, 12e-9]);
xlabel('Time (s)');
ylabel('Force (N)');
title('Comparison of Original Data vs. Imperfect Theoretical Fit');
legend('show', 'Location', 'northeast');
grid off;
hold off;

% --- 6. Save the Results ---
final_data = table(time, original_filtered_force, imperfect_force_data, imperfect_theoretical_fit, 'VariableNames', {'Time_s', 'Original_Experimental_Data', 'Artificially_Noisy_Data', 'Imperfect_Theoretical_Fit'});
writetable(final_data, 'imperfect_fit_comparison.xlsx');

disp('New Excel file "imperfect_fit_comparison.xlsx" has been created.');