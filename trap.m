%% Trapezoidal Method
% computes the step y_k+1 = y_k + h(f(t_k+1, y_k+1) + f(t_k, y_k))/2 using fixed point iteration
%
% Inputs:
% t_0: starting value for `t`
% y_0: vector of initial starting conditions for `y`
% h: step size
% t_max: final value of `t`
% f: ODE as a system of first-order equations
% J: Jacobian of `f`
%
% Outputs:
% ys: vector of the predicted values of the vector y with the ODE 
% ts: vector of time steps
function [ys, ts] = trap(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-8;

    for k = 1:length(ts)-1
        % use the previous value of y as the guess for the next one
        ys(:, k+1) = ys(:, k);
        ys(:, k+1) = fixed_point(ts(k), ys(:, k), ts(k+1), ys(:, k+1), h, tol, f);
    end
end

function [y_next] = fixed_point(t_k, y_k, t, y0, h, tol, f)
    max_iters = 100;
    iters = 0;
    % calculate the result of doing a forward euler step
    fem_value = h*f(t_k, y_k);
    y_next = fem_value;

    while norm(y_next - (y0 + h*(f(t, y_next) + fem_value)), inf) > tol && iters < max_iters
        y_next = y0 + h*(f(t, y_next) + fem_value);
        iters = iters + 1;
    end
end
