function results = fit_fractional_kv(t, y, Nu, theta_deg, p0, lb, ub, options)
% FIT_FRACTIONAL_KV Fits fractional Kelvin-Voigt constitutive model to AFM force-time data.
%
% Syntax:
%   results = fit_fractional_kv(t, y, Nu, theta_deg, p0, lb, ub, options)
%
% Inputs:
%   t         - Time vector (s)
%   y         - Measured force vector (N)
%   Nu        - Indentation velocity (m/s) [default: 2.0e-6]
%   theta_deg - Half-opening angle of conical indenter in degrees [default: 20]
%   p0        - Initial parameter guess [E0, tau, alpha] [default: [5000, 0.5, 0.05]]
%   lb        - Lower parameter bounds [default: [0, 1e-8, 0]]
%   ub        - Upper parameter bounds [default: [inf, inf, 1.0]]
%   options   - Optimset options struct [optional]
%
% Outputs:
%   results   - Struct containing optimal parameters, metrics, and fitted curve.

    if nargin < 3 || isempty(Nu), Nu = 2.0e-6; end
    if nargin < 4 || isempty(theta_deg), theta_deg = 20.0; end
    if nargin < 5 || isempty(p0), p0 = [5000, 0.5, 0.05]; end
    if nargin < 6 || isempty(lb), lb = [0, 1e-8, 0]; end
    if nargin < 7 || isempty(ub), ub = [inf, inf, 1.0]; end
    if nargin < 8 || isempty(options)
        options = optimset('TolX', 1e-6, 'TolFun', 1e-20, 'Display', 'none', 'MaxFunEvals', 3000);
    end

    t = reshape(t, [], 1);
    y = reshape(y, [], 1);
    n = length(t);
    Tr = t(end) / 2;
    theta_rad = deg2rad(theta_deg);
    KG = 1 / (2 * tan(theta_rad));

    % Objective function with exterior quadratic penalty for parameter bounds
    obj_fun = @(p) penalty_sse(p, t, y, Nu, Tr, KG, lb, ub);

    % Optimization
    [p_opt, fval, exitflag, output] = fminsearch(obj_fun, p0, options);
    p_opt = max(min(p_opt, ub), lb);

    % Compute model prediction and goodness-of-fit statistics
    y_pred = fractional_kv_force(t, p_opt, Nu, Tr, KG);
    SSres = sum((y - y_pred).^2);
    SStot = sum((y - mean(y)).^2);
    if SStot < 1e-30, SStot = 1; end
    
    R2 = 1 - (SSres / SStot);
    k = 3; % Number of fitted parameters
    adj_R2 = 1 - ((1 - R2) * (n - 1)) / (n - k - 1);
    rmse = sqrt(SSres / n);

    % Package structured output
    results.E0 = p_opt(1);
    results.tau = p_opt(2);
    results.alpha = p_opt(3);
    results.p_opt = p_opt;
    results.R_squared = R2;
    results.adjusted_R_squared = adj_R2;
    results.SSE = SSres;
    results.RMSE = rmse;
    results.y_pred = y_pred;
    results.time = t;
    results.force_exp = y;
    results.iterations = output.iterations;
    results.funcCount = output.funcCount;
end

function cost = penalty_sse(p, t, y, Nu, Tr, KG, lb, ub)
    y_pred = fractional_kv_force(t, p, Nu, Tr, KG);
    sse = sum((y - y_pred).^2);
    
    % Penalty for bound violations
    viol_lb = max(0, lb - p);
    viol_ub = max(0, p - ub);
    penalty = 1e6 * sum((viol_lb + viol_ub).^2);
    cost = sse + penalty;
end
