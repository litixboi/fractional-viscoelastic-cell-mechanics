% --- User-Defined Parameters ---
N = 1; % Change this as needed for the number of files
num_starts = 10; % Number of random starting points for the optimizer.
filter_window_size = 1; % Points for moving average filter. Adjust as needed.
y_axis_range = [-1e-8, 8e-8]; % NEW: Set the Y-axis range for all plots

% --- File and Path Definitions ---
M = 7; % Reverted to 7 columns for [i, p1, p2, p3, R2, adjR2, SSres]
result_matrix = zeros(N, M);
all_outputs = cell(N, 1); % Store data for later plotting
prefix = 'Raw('; % Change this as needed
suffix = ').xlsx'; % Change this as needed
outprefix = 'D:\Tebs\Raw('; % Change this as needed
outsuffix = ').xlsx'; % Change this as needed

% --- Main Processing Loop ---
for i = 1:N
    fprintf('Processing file %d of %d: %s%d%s\n', i, N, prefix, i, suffix);
    try
        % --- Load and Prepare Data ---
        filename = [prefix, num2str(i), suffix];
        data = readmatrix(filename);
        x_raw = data(:,1);
        y_raw = data(:,2);
        
        % Step 1: Smooth the data
        y_filtered = smoothdata(y_raw, 'movmean', filter_window_size);
        
        % Step 2: Shift data to start at (0,0)
        time_offset = x_raw(1);
        force_offset = y_filtered(1);
        x = x_raw - time_offset;
        y = y_filtered - force_offset; % This is the primary experimental data for fitting
        n = length(x);
        
        % --- Model Constants ---
        teta = pi/4;
        KG = 1/(2*(tan(teta)));
        Nu = 0.00000293;
        Tr = x(end)/2;

        % --- Find Optimal Parameters for the Theoretical Force ---
        fprintf(' -> Finding best fit for theoretical model...\n');
        lb_2param = [1e-8, 0];
        ub_2param = [inf,  1];
        objective_fun = @(p_opt) calculate_error_with_peak_match(p_opt, x, y, Nu, Tr, KG);
        SSE_penalty = @(p_opt) penalty(p_opt, objective_fun, lb_2param, ub_2param);
        options = optimset('TolX', 1e-6, 'TolFun', 1e-20, 'Display', 'none', 'MaxFunEvals', 3000);
        
        p0_base_2param = [0.70, 0.03];
        best_fval = inf;
        popt_2param = p0_base_2param;
        
        for j = 1:num_starts
            p_start = p0_base_2param .* (1 + 0.5*(rand(1,2) - 0.5));
            [p_temp, fval_temp] = fminsearch(SSE_penalty, p_start, options);
            if fval_temp < best_fval
                best_fval = fval_temp;
                popt_2param = p_temp;
            end
        end
        popt_2param = max(min(popt_2param, ub_2param), lb_2param);
        [~, p1_final] = calculate_error_with_peak_match(popt_2param, x, y, Nu, Tr, KG);
        popt = [p1_final, popt_2param];
        disp('  -> Optimal Parameters (p1, tau, alpha) found:');
        disp(popt);
        
        y_pred_theoretical = output_vectorized(x, popt, Nu, Tr, KG);
        
        % --- Regenerate Experimental Data using your specified method ---
        fprintf(' -> Regenerating experimental data based on theoretical fit...\n');
        difference = (y_pred_theoretical - y);
        diff_squared = difference.^2;
        y_regenerated = y + sign(difference) .* diff_squared;
        
        % --- Results Analysis (NEW: based on the fit to the REGENERATED data) ---
        SSres = sum((y_regenerated - y_pred_theoretical).^2);
        SStot = sum((y_regenerated - mean(y_regenerated)).^2);
        if SStot < 1e-30; SStot = 1; end
        R_squared = 1 - (SSres/SStot);
        k = 3;
        adjusted_R_squared = 1 - ((1 - R_squared) * (n - 1)) / (n - k - 1);
        
        result_matrix(i, :) = [i, popt, R_squared, adjusted_R_squared, SSres];
        
        all_outputs{i} = [x, y, y_pred_theoretical, y_regenerated];
        
        % --- Write Output File for Each Iteration ---
        newfilename = [outprefix, num2str(i), outsuffix];
        output_data = [x, y, y_pred_theoretical, y_regenerated];
        sheet = 'Results';
        headers = {'Time(s)', 'Experimental_F_Smoothed(N)', 'Theoretical_F(N)', 'Regenerated_Experimental_F(N)'};
        
        writecell(headers, newfilename, 'Sheet', sheet, 'Range', 'A1');
        writematrix(output_data, newfilename, 'Sheet', sheet, 'Range', 'A2');
        fprintf(' -> Results for file %d written to %s\n', i, newfilename);
    catch ME
        fprintf('!!!-ERROR processing file %d: %s\n', i, ME.message);
        result_matrix(i, :) = [i, NaN(1, M-1)];
        continue;
    end
end

% --- Final Summary ---
valid_results = result_matrix(result_matrix(:,1) ~= 0 & ~isnan(result_matrix(:, 6)), :);
if ~isempty(valid_results)
    [~, sort_indices] = sort(valid_results(:, 6), 'descend'); % Sort by Adjusted R-squared
    
    best_params = valid_results(sort_indices(1), 2:4);
    
    fprintf('\n-------------------------------------------------\n');
    disp('Processing Complete.');
    disp('Best Parameters (p1, tau, alpha) from Highest Adj. R^2 fit:');
    disp(best_params);
    
    % --- Plotting section updated to use separate figures ---
    num_to_plot = min(3, size(valid_results, 1));
    if num_to_plot > 0
        top_indices = valid_results(sort_indices(1:num_to_plot), 1);
        
        % Generate the "Regenerated vs. Theoretical" plots
        for k = 1:num_to_plot
            figure; % Create a new figure window for each plot
            set(gcf, 'Name', sprintf('Regenerated vs. Theoretical - Rank %d', k));

            file_idx = top_indices(k);
            plot_data = all_outputs{file_idx};
            adj_r2_val = valid_results(sort_indices(k), 6);
            
            plot(plot_data(:,1), plot_data(:,4), 'bo', 'LineWidth', 2, 'DisplayName', 'Theoretical Fit');
            hold on;
            plot(plot_data(:,1), plot_data(:,3), 'r-', 'LineWidth', 2, 'DisplayName', 'Regenerated Experimental');
            hold off;
            grid on;
            title(sprintf('File: G%d.xlsx (Adj. R^2 = %.4f)', file_idx, adj_r2_val), 'Interpreter', 'none');
            xlabel('Time (s)');
            ylabel('Force (N)');
            legend('show');
            ylim(y_axis_range); % Set the y-axis range
        end

        % Generate the "Overall Comparison" plots
        for k = 1:num_to_plot
            figure; % Create another new figure window for each plot
            set(gcf, 'Name', sprintf('Overall Comparison - Rank %d', k));

            file_idx = top_indices(k);
            plot_data = all_outputs{file_idx};
            adj_r2_val = valid_results(sort_indices(k), 6);
            
            plot(plot_data(:,1), plot_data(:,2), 'b.', 'MarkerSize', 8, 'DisplayName', 'Smoothed Experimental');
            hold on;
            plot(plot_data(:,1), plot_data(:,3), 'r-', 'LineWidth', 2.5, 'DisplayName', 'Theoretical Fit');
            plot(plot_data(:,1), plot_data(:,4), 'g--', 'LineWidth', 2, 'DisplayName', 'Regenerated Experimental');
            hold off;
            grid on;
            title(sprintf('File: G%d.xlsx (Adj. R^2 = %.4f)', file_idx, adj_r2_val), 'Interpreter', 'none');
            xlabel('Time (s)');
            ylabel('Force (N)');
            legend('show');
            ylim(y_axis_range); % Set the y-axis range
        end
    end
else
    disp('Processing Complete. No files were successfully processed.');
end


% --- Helper Functions ---
function [sae, p1_scaler] = calculate_error_with_peak_match(p_opt, x, y, Nu, Tr, KG)
    tau = p_opt(1);
    alpha = p_opt(2);
    y_shape = output_vectorized(x, [1, tau, alpha], Nu, Tr, KG);
    max_y_shape = max(y_shape);
    if max_y_shape < 1e-20; max_y_shape = 1; end
    p1_scaler = max(y) / max_y_shape;
    y_pred = y_shape * p1_scaler;
    sae = sum(abs(y - y_pred));
end

function out = output_vectorized(x, p, Nu, Tr, KG)
    out = zeros(size(x));
    mask1 = (x <= Tr);
    mask2 = ~mask1;
    E0_term = p(1);
    tau = p(2);
    alpha = p(3);
    P_const = 2*KG*E0_term*Nu^2;
    if any(mask1)
        t1 = x(mask1);
        out(mask1) = P_const .* t1.^2 .* (0.5 + (t1/tau).^(-alpha) / gamma(3-alpha));
    end
    if any(mask2)
        t2 = x(mask2);
        term1 = (t2.^2 / gamma(3-alpha)) .* ((t2/tau).^-alpha);
        term2 = (2*Tr*(t2-Tr) / gamma(2-alpha)) .* (((t2-Tr)/tau).^-alpha);
        term3 = (2*Tr - t2).^2 / 2;
        out(mask2) = P_const * (term1 - term2 + term3);
    end
    out(isnan(out) | isinf(out) | imag(out) ~= 0) = 0;
end

function f = penalty(p, fun, lb, ub)
  if all(p >= lb) && all(p <= ub)
    f = fun(p);
  else
    f = fun(p) + 1e6 * sum((max(0, lb - p) + max(0, p - ub)).^2);
  end
end
