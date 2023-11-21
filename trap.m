function [ys, ts] = trap(t_0, y_0, h, t_max, func, maxiter, tol)
    ts = t_0:h:t_max;
    ys = [y_0; zeros(t_max/h, 1)];
    for k = 1:length(ts)-1
        i = 0;
        while abs(ys(k) + h*(func(ys(k), ts(k)) + func(ts(k+1), ys(k+1))/2) - ys(k+1)) > tol && i < maxiter
            ys(k+1) = ys(k) + h*((func(ys(k), ts(k)) + func(ts(k+1), ys(k+1))/2));
            i = i + 1;
        end
    end
end
