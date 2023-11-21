function [ys, ts] = FEM(t_0, y_0, h, t_max, func)
    ts = t_0:h:t_max;
    ys = [y_0; zeros(t_max/h, 1)];
    
    for k = 1:length(ts)-1
        ys(k+1) = ys(k) + h*func(ts(k), ys(k));
    end
end
