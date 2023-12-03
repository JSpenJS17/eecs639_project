function [ys, ts] = BEM(t_0, y_0, h, t_max, f, J)
    n = length(y_0);
    ts = t_0:h:t_max;
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    tol = 1e-6;

    for k = 1:n
        % use the previos value of y as the guess for the next one
        ys(:, k+1) = newtons_method(ts(k), ys(:, k), ys(:, k), tol, f, J);
    end
end

function [y_next] = newtons_method(t, y_k, y_guess, tol, f, J)
    max_iters = 100;
    iters = 0;
    y_next = y_guess;

    while norm(f(t, y_next), inf) > tol && iters < max_iters
        J_tyk = J(t, y_next);
        % replace this with LUP decomp when i actually get it working
        s_k = J_tyk \ -f(t, y_next);

        y_next = y_k + s_k;
        iters = iters + 1;
    end
end
