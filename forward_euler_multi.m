function [t, y] = eulers_method_multi(h, t_0, t_max, y_0, u)
    n = length(y_0);
    t = linspace(t_0, t_max, (t_max - t_0) / h);
    y = [y_0 zeros(length(y_0), length(t) - 1)];

    for k = 2:length(t)
        for i = 1:n - 1
            y(i, k) = y(i, k-1) + h * y(i+1, k-1);
        end

        y(n, k) = y(n, k-1) + h * u(t(k), y(:, k));
    end
end