function [ys, ts] = FEM(t_0, y_0, h, t_max, f)
    % detailed preamble comment goes here
    ts = (t_0:h:t_max)';
    ys = [y_0;zeros(length(ts)-1, length(y_0))];

    for k = 1:length(ts)-1
        ys(k+1, :) = ys(k, :) + h*f(ts(k, :), ys(k, :));

    end
end
