%% TEST VIBRATING SPRING
% exact solution: y = .2*cos(8*t)

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [.2; 0];
t_max = 5;

% define exact solution using symbolic notation for comparison later
syms t
y = .2*cos(8*t);

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @f, y, "Vibrating Spring", "Meters")

%% TEST ELECTRIC CIRCUIT
% exact solution: (4/697)*((exp(-20)*t/3)*(-63*cos(15*t)-116*sin(15*t)) + (21*cos(10*t) + 16*sin(10*t)))

% define starting conditions from project description
h = .025;
t_0 = 0;
y_0 = [0; 0];
t_max = 5;

% define exact solution using symbolic notation for comparison later
syms t
y = (4/697)*((exp(-20)*t/3)*(-63*cos(15*t)-116*sin(15*t)) + (21*cos(10*t) + 16*sin(10*t)));

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @u, y, "Electric Circuit", "Coulombs")

%% FUNCTION DEFINITIONS
% define f for vibrating spring
function [z] = f(t, y)
    % z = [u2; -kx/m], from project description
    z = [y(2); -64*y(1)];
end

% define u for electric circuit
function [z] = u(t, y)
    % z = [u2; E(t) - RQ' - (1/C)Q
    R = 40;
    C = 16*10^-4;
    z = [y(2); 100*cos(10*t) - R*y(2) - (1/C)*y(1)];
end

% define helper function to run solvers and plot results
function [] = plot_test(t_0, y_0, h, t_max, f, exact_soln, test_name, y_label)
    figure;
    
    disp("Results in order: FEM, BEM, ITM, RK4, ABM4");

    % run FEM
    tic
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,1);
    plot(ts, ys(1, :), 'r.-');
    xlabel("Time (t)");
    ylabel(y_label);
    title('FEM of ' + test_name);
    
    % run BEM,  disabled until root finding is finished
    tic
    [ys, ts] = BEM(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,2);
    plot(ts, ys(1, :), 'g.-');
    xlabel("Time (t)");
    ylabel(y_label);
    title('BEM of ' + test_name);
  
    % run trap on electric circuit, disabled until fixed
    tic
    [ys, ts] = trap(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,3);
    plot(ts, ys(1, :), 'b.-');
    xlabel("Time (t)");
    ylabel(y_label);
    title('ITM of ' + test_name);
    
    % run RK4
    tic
    [ys, ts] = RK4(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,4);
    plot(ts, ys(1, :), 'm.-');
    xlabel("Time (t)");
    ylabel(y_label);
    title('RK4 of ' + test_name);
    
    % run AB4 and AM4
    tic
    [ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,5);
    plot(ts, ys(1, :), 'k.-');
    xlabel("Time (t)");
    ylabel(y_label);
    title('ABM4 of ' + test_name);
    fprintf("\n");

    % plot the exact solution for comparison
    subplot(3,2,6);
    fplot(exact_soln, [t_0 t_max], 'c.-');
    ylim([-.2 .2]);
    title('Exact Solution of ' + test_name);
end
