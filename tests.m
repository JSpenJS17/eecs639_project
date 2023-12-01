hold on;
clf;

%% y'' = -y with y(0) = 0, y'(0) = 1
%
% u1 = y
% u2 = y'
%
% u1' = u2
% u2' = -u1

h = .01;
t_0 = 0;
y_0 = [0; 1];
t_max = pi*2;

f = @(t, y) [y(2); -y(1)];

%% Testing FEM
[ys, ts] = FEM(t_0, y_0, h, t_max, f);
% plot(ts, ys(1, :), 'black');

%% Testing RK4
[ys, ts] = RK4(t_0, y_0, h, t_max, f);
% plot(ts, ys(1, :), 'black');

%% Testing Predictor-Corrector
[ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f);
plot(ts, ys(1, :), 'black');
