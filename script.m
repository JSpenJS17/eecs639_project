%% TEST VIBRATING SPRING
% exact solution: y = .2*cos(8*t)

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [.2; 0];
t_max = 5;
syms t
y = .2*cos(8*t);

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @f, "Vibrating Spring", true, y)

%% TEST ELECTRIC CIRCUIT
% exact solution: (4/697)*((exp(-20)*t/3)*(-63*cos(15*t)-116*sin(15*t)) + (21*cos(10*t) + 16*sin(10*t)))

% define starting conditions from project description
h = .01;
t_0 = 0;
y_0 = [0; 0];
t_max = 5;
syms t
y = (4/697)*((exp(-20)*t/3)*(-63*cos(15*t)-116*sin(15*t)) + (21*cos(10*t) + 16*sin(10*t)));

%test the different IVP methods
plot_test(t_0, y_0, h, t_max, @u, "Electric Circuit", true, y)


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
function [] = plot_test(t_0, y_0, h, t_max, f, test_name, plot_exact, exact_soln)
    figure;

    % run FEM
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    subplot(3,2,1);
    plot(ts, ys(1, :), 'r.-');
    title('FEM of ' + test_name)
    
    % (HARLAN) %
    % run BEM,  disabled until root finding is finished
    %[ys, ts] = BEM(t_0, y_0, h, t_max, f);
    %subplot(3,2,2);
    %plot(ts, ys(1, :), 'g.-');
    %title('BEM of ' + test_name)
  
    % (PIERCE) %
    % run trap on electric circuit, disabled until fixed
    %[ys, ts] = trap(t_0, y_0, h, t_max, f);
    %subplot(3,2,3);
    %plot(ts, ys(1, :), 'b.-');
    %title('Trap of ' + test_name)
    
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
    
    %plot exact solution (optional)
    if plot_exact
        subplot(3,2,6);
        fplot(exact_soln, [t_0 t_max], 'c.-')
        title('Exact Solution of ' + test_name)
    end

end
