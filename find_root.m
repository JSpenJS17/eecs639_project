% Use Newton's Method to find a root of the function `f`
% with first derivative `fp`, starting from the point `x`
function [x] = find_root(f, fp, x)
    tol = 1e-8;
    max_iters = 100;
    k = 0;

    while abs(f(x)) > tol && k < max_iters
        x = x - (f(x) / fp(x));
        k = k + 1;
    end
end
