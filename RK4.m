% Runge-Kutta 4th order method
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
