% Predictor-corrector method using AB4 as a predictor
% and AM4 as a corrector.
function [ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    ys(:, 2) = RK4_step(ts(1), ys(:, 1), h, f);
    ys(:, 3) = RK4_step(ts(2), ys(:, 2), h, f);
    ys(:, 4) = RK4_step(ts(3), ys(:, 3), h, f);

    for k = 4:length(ts)-1
        ys(:, k+1) = AB4_predictor(k, ts, ys, h, f);
        ys(:, k+1) = AM4_corrector(k, ts, ys, h, f);
    end
end

function [y] = AB4_predictor(k, ts, ys, h, f)
    y = ys(:, k) + (h / 24) * (55*f(ts(k), ys(:, k)) - 59*f(ts(k-1), ys(:, k-1)) + 37*f(ts(k-2), ys(:, k-2)) - 9*f(ts(k-3), ys(:, k-3)));
end

function [y] = AM4_corrector(k, ts, ys, h, f)
    y = ys(:, k) + (h / 24) * (9*f(ts(k+1), ys(:, k+1)) + 19*f(ts(k), ys(:, k)) - 5*f(ts(k-1), ys(:, k-1)) + f(ts(k-2), ys(:, k-2)));
end

function [y] = RK4_step(t, y_0, h, f)
    k1 = f(t,       y_0);
    k2 = f(t + h/2, y_0 + k1*(h/2));
    k3 = f(t + h/2, y_0 + k2*(h/2));
    k4 = f(t + h,   y_0 + k3*h);
    y = y_0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end
