% Define the number of files to process
N = 1; % Change this as needed
M = 7;
result_matrix = zeros(N, M);
% Define the common part of the file names
prefix = 'Raw('; % Change this as needed
suffix = ').xlsx'; % Change this as needed
% Define the output file name format
outprefix = 'D:\Tebs\Raw('; % Change this as needed
outsuffix = ').xlsx'; % Change this as needed
% Loop over the files
for i = 1:N
    fprintf('Processing file %d of %d: %s%d%s\n', i, N, prefix, i, suffix);

    try
        % Construct the input file name
        filename = [prefix, num2str(i), suffix];
        
        % Read the excel file
        data = readmatrix(filename);
        
        % --- Adjust Data to Start from Zero ---
        x_raw = data(:,1);
        y_raw = data(:,2);
        time_offset = x_raw(1);
        force_offset = y_raw(1);
        x = x_raw - time_offset;
        y = y_raw - force_offset;
        
        % --- Y-Data Scaling for Numerical Stability ---
        y_scale_factor = max(y);
        if y_scale_factor < 1e-20 % Avoid division by zero if force is zero
            y_scale_factor = 1;
        end
        y_scaled = y / y_scale_factor;

        % --- Calculations based on shifted data ---
        n = length(x);
        teta = pi/4;
        KG = 1/(2*(tan(teta)));
        Nu = 0.00000293;
        Tr = x(end)/2;

        % --- Optimization using a more ACCURATE fminsearch strategy ---
        lb = [0, 1e-8, 0]; % Lower bounds
        ub = [inf, inf, 1]; % Upper bounds
        
        % The objective function now works with SCALED data
        sse_fun = @(p) sum((y_scaled - (output_vectorized(x, p, Nu, Tr, KG) / y_scale_factor)).^2);
        
        % Define the SSE with the penalty function from your original code
        SSE_penalty = @(p) penalty(p, sse_fun, lb, ub);
        
        % Set options for fminsearch
        options = optimset('TolX', 1e-6, 'TolFun', 1e-20, 'Display', 'none', 'MaxFunEvals', 2000);

        % --- NEW: Multi-start strategy for improved accuracy ---
        % Run the optimizer multiple times from different starting points
        % to avoid getting stuck in a local minimum.
        num_starts = 5;
        p0_base = [15000, 0.70, 0.03]; % A neutral base starting point
        
        best_fval = inf; % Best error found so far, start at infinity
        popt = p0_base; % Best parameters found so far
        
        fprintf(' -> Running fminsearch with %d starting points...\n', num_starts);
        for j = 1:num_starts
            % Create a slightly randomized starting point for each run
            p_start = p0_base .* (1 + 0.5*(rand(1,3) - 0.5)); % +/- 25% random variation
            
            [p_temp, fval_temp] = fminsearch(SSE_penalty, p_start, options);
            
            % If this run found a better solution, keep it
            if fval_temp < best_fval
                best_fval = fval_temp;
                popt = p_temp;
            end
        end
        
        % Project final parameters to ensure they are within bounds
        popt = max(min(popt, ub), lb);
        
        % Display the optimal parameters
        disp('Optimal Parameters (popt) found by robust fminsearch:');
        disp(popt);
        
        % Calculate predicted values using the optimal parameters
        y_pred = output_vectorized(x, popt, Nu, Tr, KG);
        
        % --- Results Analysis ---
        SSres = sum((y - y_pred).^2);
        SStot = sum((y - mean(y)).^2);
        if SStot < 1e-30; SStot = 1; end % Avoid division by zero
        R_squared = 1 - (SSres/SStot);
        
        k = 3; % Number of predictor variables
        adjusted_R_squared = 1 - ((1 - R_squared) * (n - 1)) / (n - k - 1);
         
        % Store the results in the result matrix
        result_matrix(i, :) = [i, popt, R_squared, adjusted_R_squared, SSres];
        
        % --- Write Output File for Each Iteration ---
        newfilename = [outprefix, num2str(i), outsuffix];
        output_data = [x, y, y_pred];
        sheet = 'Results';
        headers = {'Time_shifted(s)', 'Experimental_F_shifted(N)', 'Theoretical_F(N)'};
        
        writecell(headers, newfilename, 'Sheet', sheet, 'Range', 'A1');
        writematrix(output_data, newfilename, 'Sheet', sheet, 'Range', 'A2');
        
        fprintf(' -> Results for file %d written to %s\n', i, newfilename);

    catch ME
        % --- NEW: More informative error catching ---
        fprintf('!!!-ERROR processing file %d: %s\n', i, ME.message);
        result_matrix(i, :) = [i, NaN(1, 6)]; % Mark row as invalid
        continue; % Skip to the next file
    end
end

% --- Final Summary (outside the loop) ---
% Filter out failed runs before calculating the best result
valid_results = result_matrix(~isnan(result_matrix(:, 6)), :);
if ~isempty(valid_results)
    [~, idx] = max(valid_results(:, 6)); % Find max adj R-squared from valid results
    best_values_R = valid_results(idx, 2:4);
    fprintf('\n-------------------------------------------------\n');
    disp('Processing Complete.');
    disp('Best set of values (Highest Adjusted R^2):');
    disp(best_values_R);

    [~, idx] = min(valid_results(:, 7)); % Find min SSres from valid results
    best_values_min = valid_results(idx, 2:4);
    disp('Best set of values (Lowest SSres):');
    disp(best_values_min);
else
    fprintf('\n-------------------------------------------------\n');
    disp('Processing Complete. No files were successfully processed.');
end


% --- Helper Functions ---

function out = output_vectorized(x, p, Nu, Tr, KG)
    % This is the fast, vectorized version of your original output function
    out = zeros(size(x));
    
    mask1 = (x <= Tr);
    mask2 = ~mask1;
    
    P_const = 2*KG*p(1)*Nu^2;
    tau = p(2);
    alpha = p(3);

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
  % This penalty function is from your original code.
  % It adds a large penalty if parameters go out of bounds.
  if all(p >= lb) && all(p <= ub)
    f = fun(p);
  else
    f = fun(p) + 1e6 * sum((max(0, lb - p) + max(0, p - ub)).^2);
  end
end