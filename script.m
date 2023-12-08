%% TEST VIBRATING SPRING
% exact solution: y = .2*cos(8*t)

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [.2; 0];
t_max = 5;

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @f, "Vibrating Spring")

%% TEST ELECTRIC CIRCUIT
% exact solution: (4/697)*((exp(-20)*t/3)*(-63*cos(15*t)-116*sin(15*t)) + (21*cos(10*t) + 16*sin(10*t)))

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [0; 0];
t_max = 5;

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @u, "Electric Circuit")

%% DO OREGONATOR

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [0; 0.23; 0.0];
t_max = 25;

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @o, "Oregonator")

%% FUNCTION DEFINITIONS
% define f for vibrating spring
function [z] = f(t, y)
    % z = [u2; -kx/m], from project description
    z = [y(2); -64*y(1)];
end

% define Jacobian of f for BEM
function [J] = Jf(t, y)
    J = [
        0 1;
        -64 0
    ];
end

% define u for electric circuit
function [z] = u(t, y)
    % z = [u2; E(t) - RQ' - (1/C)Q
    R = 40;
    C = 16*10^-4;
    z = [y(2); 100*cos(10*t) - R*y(2) - (1/C)*y(1)];
end

% define Jacobian of u for BEM
function [J] = Ju(t, y)
    R = 40;
    C = 16*10^-4;
    J = [
        0 1;
        -1/C -R
    ];
end

% define o for oregonator
function [u] = o(t, y)
    % reaction parameters
    A = 0.06;
    B = 0.02;
    f = 1.0;
    k3 = 1.0;
    k2 = 2.4e5;
    k5 = 1.28;
    k4 = 3.0e3;
    k0 = 33.6;

    % X, Y, and Z
    X = y(1);
    Y = y(2);
    Z = y(3);

    u = [
            ( k3*A*Y - k2*X*Y + k5*A*X);
            (-k3*A*Y - k2*X*Y + .5*f*k0*B*Z );
            ( 2*k5*A*X - k0*B*Z )
        ];
end

% define Jacobian of o for BEM
function [J] = Jo(t, y)
    % reaction parameters
    a = 0.06;
    b = 0.02;

    kc = 1.0;
    k2 = 2.4E+06;
    k3 = 1.28;
    k4 = 3.0E+03;
    k5 = 33.6;

    %calculated scaling
    eta1 = kc * b / k5 / a;
    eta2 = 2.0 * kc * k4 * b / k2 / k5 / a;
    q = 2.0 * k3 * k4 / k2 / k5;
    f = 1.0;

    J = [
        ((-y(2) + 1 - 2.*y(1))/eta1) (q - y(1))/eta1    0;
        -y(2)/eta2                   (-q - y(1))/eta2   f/eta2;
        1                            0                  -1
    ];
end

% define helper function to run solvers and plot results
function [] = plot_test(t_0, y_0, h, t_max, f, test_name)
    figure;

    % run FEM
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    subplot(3,2,1);
    plot(ts, ys(1, :), 'r.-');
    title('FEM of ' + test_name)
    
    % run BEM,  disabled until root finding is finished
    [ys, ts] = BEM(t_0, y_0, h, t_max, f);
    subplot(3,2,2);
    plot(ts, ys(1, :), 'g.-');
    title('BEM of ' + test_name)
  
    % run trap on electric circuit, disabled until fixed
    [ys, ts] = trap(t_0, y_0, h, t_max, f);
    subplot(3,2,3);
    plot(ts, ys(1, :), 'b.-');
    title('Trap of ' + test_name)
    
    % run RK4
    [ys, ts] = RK4(t_0, y_0, h, t_max, f);
    subplot(3,2,4);
    plot(ts, ys(1, :), 'm.-');
    title('RK4 of ' + test_name)
    
    % run AB4 and AM4
    [ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f);
    subplot(3,2,5);
    plot(ts, ys(1, :), 'k.-');
    title('PC of ' + test_name)

end
