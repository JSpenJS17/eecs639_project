function [ys, ts] = trap(t_0, y_0, h, t_max, f, J)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-6;

    for k = 1:length(ts)-1
        % use the previous value of y as the guess for the next one
        ys(:, k+1) = ys(:, k);
        % ys(:, k+1) = newtons_method(ts(k+1), ys(:, k+1), h, tol, f, J);
        ys(:, k+1) = fixed_point(ts(k), ys(:, k), ts(k+1), ys(:, k+1), h, tol, f);
    end
end

function [y_next] = fixed_point(t_k, y_k, t, y0, h, tol, f)
    max_iters = 100;
    iters = 0;
    fem_value = h*f(t_k, y_k);
    y_next = fem_value;

    while norm(y_next - (y0 + h*f(t, y_next) + fem_value)/2, inf) > tol && iters < max_iters
        y_next = y0 + h*(f(t, y_next) + fem_value)/2;
        iters = iters + 1;
    end
end
