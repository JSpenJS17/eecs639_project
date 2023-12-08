%% Forward Euler's Method
% solves the step y_k+1 = y_k + hf(t_k, y_k)
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
function [ys, ts] = FEM(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = (t_0:h:t_max)';
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;

    for k = 1:length(ts)-1
        ys(:, k+1) = ys(:, k) + h*f(ts(k), ys(:, k));
    end
end