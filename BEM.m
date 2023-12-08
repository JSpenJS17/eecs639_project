%% Backward Euler's Method
% solves the step y_k+1 = y_k + hf(t_k+1, y_k+1) using both
% Newton's method and fixed-point iteration when necessary.
%
% Inputs:
% t_0: starting value for `t`
% y_0: vector of initial starting conditions for `y`
% h: step size
% t_max: final value of `t`
% f: ODE as a system of first-order equations
%
% Outputs:
% ys: vector of the predicted values of the vector y with the ODE 
% ts: vector of time steps
function [ys, ts] = BEM(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-8;

    for k = 1:length(ts)-1
        % use the previous value of y as the guess for the next one
        ys(:, k+1) = ys(:, k);
        ys(:, k+1) = fixed_point(ts(k+1), ys(:, k+1), h, tol, f);
    end
end

function [y_next] = fixed_point(t, y0, h, tol, f)
    max_iters = 100;
    iters = 0;
    y_next = y0;

    % do fixed point iteration to find y_next
    while norm(y_next - (y0 + h*f(t, y_next)), inf) > tol && iters < max_iters
        y_next = y0 + h*f(t, y_next);
        iters = iters + 1;
    end
end
