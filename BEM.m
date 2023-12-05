function [ys, ts] = BEM(t_0, y_0, h, t_max, f, J)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-6;

    for k = 1:n
        % use the previos value of y as the guess for the next one
        ys(:, k+1) = ys(:, k);
        ys(:, k+1) = newtons_method(ts(k+1), ys(:, k+1), h, tol, f, J);
        % ys(:, k+1) = fixed_point(ts(k+1), ys(:, k+1), h, tol, f);
    end
end

function [y_next] = newtons_method(t, y0, h, tol, f, J)
    max_iters = 100;
    iters = 0;
    y_next = y0;

    while norm(y_next - (y0 + h*f(t, y_next)), inf) > tol && iters < max_iters
        [L, U, P] = lu(J(t, y_next));
        v = L\P'*-(y_next - y0 - h*f(t, y_next));
        s_k = U\v;

        y_next = y_next + s_k;
        iters = iters + 1;
    end
end

function [y_next] = fixed_point(t, y0, h, tol, f)
    max_iters = 100;
    iters = 0;
    y_next = y0;

    while norm(y_next - (y0 + h*f(t, y_next)), inf) > tol && iters < max_iters
        y_next = y0 + h*f(t, y_next);
        iters = iters + 1;
    end
end
