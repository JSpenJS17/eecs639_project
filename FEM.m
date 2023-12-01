% Forward Eulers method
function [ys, ts] = FEM(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;

    for k = 1:length(ts)-1
        ys(:, k+1) = ys(:, k) + h*f(ts(k), ys(:, k));
    end
end
