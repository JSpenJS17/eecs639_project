function [z] = f(t, y)
    % z = [u2, f(t, y)], like from class.
    % in this case, f(t, y) is just -y, meaning that the function is either
    % cos(x) or sin(x), depending on the initial conditions.
    z = [y(2) -y(1)];
end

