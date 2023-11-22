function [ys, ts] = BEM(t_0, y_0, h, t_max, func, maxiter, tol)
    ts = (t_0:h:t_max)';
    ys = [y_0;zeros(length(ts)-1, length(y_0))];
    
    for k = 1:length(ts)-1
        i = 0;
        while abs(norm(ys(k, :) + h*func(ts(k+1), ys(k+1, :) - ys(k+1, :)))) > tol && i < maxiter
            ys(k+1, :) = ys(k, :) + h*func(ts(k+1), ys(k+1, :));
            i = i + 1;
        end
    end
end
