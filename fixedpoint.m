function [x] = fixedpoint(g, x, t, tol, maxiter)
%does the fixed point iteration of g
%given guess x and tolerance tol
    %loop until we find a fixed point
    k = 0;
    while abs(g(x, t) - x) > tol && k < maxiter
        %set the next x = g(old x)
        x = g(x, t);
        k = k + 1;
    end
end
