%% TESTING FEM
h = .01;
t_0 = 0;
y_0 = 1;
t_max = 2;
[ys, ts] = FEM(t_0, y_0, h, t_max, @f);
figure;
plot(ts, ys, 'k.-');
syms x
g = exp(x);
figure;
fplot(g, [t_0 t_max])

%% TESTING BEM
h = .01;
t_0 = 0;
y_0 = 1;
t_max = 2;
tol = 10^-8;
maxiter = 50;
[ys, ts] = BEM(t_0, y_0, h, t_max, @f, maxiter, tol);
figure;
plot(ts, ys, 'k.-');
syms x
g = exp(x);
figure;
fplot(g, [t_0 t_max])

%% TESTING trap
h = .01;
t_0 = 0;
y_0 = 1;
t_max = 2;
tol = 10^-8;
maxiter = 50;
[ys, ts] = trap(t_0, y_0, h, t_max, @f, maxiter, tol);
figure;
plot(ts, ys, 'k.-');
syms x
g = exp(x);
figure;
fplot(g, [t_0 t_max])


% TESTING forward_euler_multi
h = 0.025;
t_0 = 0;
t_max = 3.5;
y_0 = [0; 1];

[t, y] = forward_euler_multi(h, t_0, t_max, y_0, @u);
figure;
plot(t, y(2, :));

% TESTING backward_euler_multi
h = 0.025;
t_0 = 0;
t_max = 3.5;
y_0 = [0; 1];

[t, y] = backward_euler_multi(h, t_0, t_max, y_0, @u);
figure;
plot(t, y(2, :));

% y'' = -y
function [z] = u(t, y)
    z = -y(1);
end