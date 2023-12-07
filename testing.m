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
hold on
plot(ts, ys(:,1), 'r');
plot(ts, ys(:,2), 'g');
plot(ts, ys(:,3), 'b');
hold off