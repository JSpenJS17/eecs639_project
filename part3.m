clf;
hold on;

h = .01;
t_0 = 0;
t_max = pi*2;

a = 0.06;
b = 0.02;
k0 = 1.0;
k2 = 2.4e5;
k3 = 1.28;
k4 = 3.0e3;
k5 = 33.6;
f = 1.0;

y_0 = [0.0; 0.; 0.023];

[ts, ys] = ode45(dy, [0 1], y_0);

% plot(ts, ys(:, 1), 'green');
plot(ts, ys(:,1), 'red');
% plot(ys(:,1), ys(:,2), 'green');
