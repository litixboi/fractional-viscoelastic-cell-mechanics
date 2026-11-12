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
    % Construct the input file name
    filename = [prefix, num2str(i), suffix];
    % Read the excel file and store the data in a matrix
    % Use readmatrix instead of xlsread for better performance
    data = readmatrix(filename);
    % Extract the second column (x) and assign it to a variable
    len = data(:,1);
    n = length(len);
    x = data(:,1);
    y = data(:,2);
    teta = pi/4;
    KG = 1/(2*(tan(teta)));
    Nu = 0.00000293;
    Tr = x(n)/2;
    epsilon = 0.0000001;

    % Define the lower and upper bounds for the parameters
    lb = [0, 0, 0]; % Lower bounds
    ub = [inf, inf, 1]; % Upper bounds
    % Define a function that computes the SSE with a penalty given the parameters
    SSE = @(p) penalty (p,@ (p) sum ((y - output (x, p, Nu, Tr, KG)).^2),lb,ub); % Use abs (p) to make the parameters positive

    % Set an initial guess for the parameters
    p0 = [15000, 0.70, 0.03];
    options = optimset ('TolX',1e-6,'TolFun',1e-20,'Display','iter','MaxFunEvals',2000); % increase the maximum number of function evaluations
    % Use fminsearch to find the optimal parameters without bounds
    popt = fminsearch (SSE,p0,options);


    % Display the optimal parameters
    disp(popt); % Use popt to make the parameters positive
    % Calculate the output values using the optimal parameters
    out = output(x, popt, Nu, Tr, KG); % Use popt to make the parameters positive
    % Display the standard deviation of the output values
    disp(std(out));
    % Calculate the predicted values using the optimal parameters
    y_pred = output(x, popt, Nu, Tr, KG); % Use popt to make the parameters positive
    % Calculate the R-squared and adjusted R-squared values
    SSres = sum((y - y_pred).^2);
    SStot = sum((y - mean(y)).^2);
    R_squared = 1 - (SSres/SStot);
    n = length(y);
    k = 3; % Number of predictor variables
    adjusted_R_squared = 1 - ((1 - R_squared) * (n - 1)) / (n - k - 1);
    
    % Store the results in the result matrix
    result_matrix(i, :) = [i, popt, R_squared, adjusted_R_squared, SSres];
    % Find the set of values with the highest adjusted R-squared
    [~, idx] = max(result_matrix(:, 6)); % Assuming R-squared is in the 5th column
    best_values_R = result_matrix(idx, 2:4);
    disp('Best set of values:(R^2adjusted)');
    disp(best_values_R);
    [~, idx] = min(result_matrix(:, 7)); % Assuming R-squared is in the 5th column
    best_values_min = result_matrix(idx, 2:4);
    disp('Best set of values:(min)');
    disp(best_values_min);
    % Construct the output file name
    newfilename = [outprefix, num2str(i), outsuffix];
    % Write the x value from the initial excel data to the 1st column of a new excel file
    % Write the y data from the initial excel file to the 2nd column of the new excel file
    % Write the result of the function for each variable x to the 3rd column of the new excel file
    % Specify the sheet name and the column names
    sheet = 'Results'; % Name of the sheet
    writecell({'Time(s)', 'Experimental F(N)', 'Theoretical F(N)'}, newfilename, 'Sheet', sheet); % Write the column names to the next available row, without specifying the range, using writecell instead of writematrix
    writematrix(x, newfilename, 'Sheet', sheet, 'Range', ['A2:A' num2str(length(x))]); % Write x to the 1st column, specifying the range, and transpose it
    writematrix(y, newfilename, 'Sheet', sheet, 'Range', ['B2:B' num2str(length(y))]); % Write y to the 2nd column, specifying the range, and transpose it
    writematrix(output(x, popt, Nu, Tr, KG), newfilename, 'Sheet', sheet, 'Range', ['C2:C' num2str(length(output(x, popt, Nu, Tr, KG)))]); % Write the output to the 3rd column, specifying the range, and transpose it
end


function out = output(x, p, Nu, Tr, KG) %13th eq
    % Initialize an empty vector to store the output values
    out = zeros(size(x));
    % Use parfor loop instead of while loop to parallelize the computation
    % This requires a Parallel Computing Toolbox license
    parfor i = 1:length(x)
        % Check the condition for x
        if x(i) <= Tr
            out(i) = 2*KG*p(1)*Nu^2*x(i)^2*(1/2 + (((x(i)/p(2))^(-p(3)))/gamma(3 - p(3))));
        else
            out(i) = 2*KG*p(1)*Nu^2*((x(i)^2/gamma(3 - p(3)))*((x(i)/p(2))^-p(3)) - (2*Tr*(x(i)-Tr))/gamma(2 - p(3))*(((x(i)-Tr)/p(2))^(-p(3))) + ((2*Tr - x(i))^2)/2);
        end
    end
end
function f = penalty (p,fun,lb,ub)
  % Check if the parameters are within the bounds
  if all (p >= lb) && all (p <= ub)
    % Evaluate the original objective function
    f = fun (p);
  else
    % Add a large penalty value
    f = fun (p) + 1e6;
  end
end
function p = transform (x,lb,ub)
  % Map the unbounded variables x to the bounded parameters p
  p = lb + (ub - lb)./(1 + exp (-x));
end
function p = project (p,lb,ub)
  % Project the parameters p back to the bounds
  p = max (min (p,ub),lb);
end
