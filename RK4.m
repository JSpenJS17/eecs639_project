% Runge-Kutta 4th-Order Method
% solves the step y_k+1 using a 4th order accurate method
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
function [ys, ts] = RK4(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;

    for k = 1:length(ts)-1
        k1 = f(ts(k), ys(:, k));
        k2 = f(ts(k) + h/2., ys(:, k) + k1*(h/2.));
        k3 = f(ts(k) + h/2., ys(:, k) + k2*(h/2.));
        k4 = f(ts(k) + h, ys(:, k) + k3*h);
        ys(:, k+1) = ys(:, k) + (h/6.) * (k1 + 2.*k2 + 2.*k3 + k4);
    end
end