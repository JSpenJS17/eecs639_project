function [ys, ts] = BEM(t_0, y_0, h, t_max, f, J)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-6;

    for k = 1:length(ts)-1
        % use the previous value of y as the guess for the next one
        ys(:, k+1) = ys(:, k);
        ys(:, k+1) = newtons_method(ts(k+1), ys(:, k+1), h, tol, f, J);
    end
end

function [y_next] = newtons_method(t, y0, h, tol, f, J)
    max_iters = 100;
    iters = 0;
    y_next = y0;

    error = norm(y_next - (y0 + h*f(t, y_next)), inf);

    while error > tol && iters < max_iters
        max_correction_iters = 100;
        correction_iters = 0;
        [L, U, P] = lu(J(t, y_next));
        v = L\P'*-(y_next - y0 - h*f(t, y_next));
        s_k = U\v;

        y_next = y_next + s_k;
        iters = iters + 1;
        current_error = norm(y_next - (y0 + h*f(t, y_next)), inf);
        while current_error > error && correction_iters < max_correction_iters
            % perform fixed point iter until error is reduced
            y_next = y0 + h*f(t, y_next);
            current_error = norm(y_next - (y0 + h*f(t, y_next)), inf);
            correction_iters = correction_iters + 1;
        end

        error = current_error;
    end
end
