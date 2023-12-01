function [ys, ts] = BEM(t_0, y_0, h, t_max, f)
    n = length(y_0);
    ts = (t_0:h:t_max)';
    ys = zeros(n, length(ts));
    ys(:, 1) = y_0;
    
    for k = 1:length(ts)-1
        i = 0;
        u = @(t, y) ys(:, k) + h*f(ts(k+1), ys(:, k+1));
        %this doesn't work, just needs the jacobian bit
        ys(:, k+1) = find_root_fd(u, ys(k));
    end
end

function [x] = find_root_fd(f, x)
    tol = 1e-8;
    max_iters = 100;
    k = 0;

    while norm(abs(f(x))) > tol && k < max_iters
        x = x - (f(x) / finite_difference(f, x, tol));
        k = k + 1;
    end
end

function [dx] = finite_difference(f, x, h)
    dx = (f(x + h) - f(x - h)) / (2*h);
end