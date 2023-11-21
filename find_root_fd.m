% Use Newton's Method to find a root of the function `f`
% by approximating the derivative f'(x) with the finite
% difference method, starting from the point `x`
function [x] = find_root_fd(f, x)
    tol = 1e-8;
    max_iters = 100;
    k = 0;

    while abs(f(x)) > tol && k < max_iters
        x = x - (f(x) / finite_difference(f, x, tol));
        k = k + 1;
    end
end

function [dx] = finite_difference(f, x, h)
    dx = (f(x + h) - f(x - h)) / (2*h);
end
