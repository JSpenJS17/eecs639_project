clf;
hold on;

h = .01;
t_0 = 0;
t_max = 5;

%% reaction parameters
a = 0.06;
b = 0.02;
k0 = 1.0;
k2 = 2.4e5;
k3 = 1.28;
k4 = 3.0e3;
k5 = 33.6;
f = 1.0;

y_0 = [0.0; 0.; 0.023];

dy = @(t, y) [
    ( k3*a*y(2) - k2*y(1)*y(2) + k5*a*y(1) -2*k4*y(1).^2);
    (-k3*a*y(2) - k2*y(1)*y(2) + 0.5*f*k0*b*y(3));
    (2*k5*a*y(1) - k0*b*y(3))
];

[ts, ys] = ode45(dy, [t_0 t_max], y_0);

% plot X(t) with respect to t
plot(ts, ys(:,1), 'red');

% plot concentration of Y with respect to X
% plot(ys(:,1), ys(:,2), 'green');


