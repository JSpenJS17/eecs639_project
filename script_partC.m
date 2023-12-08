% define starting time and timestep variables
h = .0000025;
t_0 = 0;
t_max = 25;

% define starting x, y, and z
y_0 = [1;1;1;];

% test all 5 plots
plot_test_o(t_0, y_0, h, t_max, @o, "BZ Model");

% DEFINE FUNCTIONS

function [dydt] = o(t, y)
    % this function works with these starting parameters, but it uses the
    % scaled version instead of the regular one given in the project

    % reaction parameters
    a = 0.06;
    b = 0.02;
    f = 1.0;

    %set k0-5 as the known constants
    k0 = 1.0;
    k2 = 2.4E+06;
    k3 = 1.28;
    k4 = 3.0E+03;
    k5 = 33.6;

    % calculate scaling stuff
    eps = k0 * b / k5 / a;
    eps2 = 2.0 * k0 * k4 * b / k2 / k5 / a;
    q = 2.0 * k3 * k4 / k2 / k5;

    % grab x, y, and z (using u, v, and w to avoid naming conflicts)
    u = y(1);
    v = y(2);
    w = y(3);
    
    % calculate new x, y, and z
    dxdt = ( q*v - u*v + u*(1.0 - u) ) / eps;
    dydt = ( -q*v - u*v + f*w ) / eps2;
    dzdt = u - w;
    
    % return them
    dydt = [ dxdt; dydt; dzdt ];

end

% define helper function to run solvers and plot results
function [] = plot_test_o(t_0, y_0, h, t_max, f, test_name)
    figure;

    disp("Results in order: FEM, BEM, ITM, RK4, ABM4");

    % run FEM
    tic
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,1);
    hold on
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');
    hold off
    xlabel("Time (t)");
    ylabel("log10(x, y, z)")
    title('FEM of ' + test_name)
    
    % run BEM
    tic
    [ys, ts] = BEM(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,2);
    hold on
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');
    hold off
    xlabel("Time (t)");
    ylabel("log10(x, y, z)")
    title('BEM of ' + test_name)
  
    % run trap
    tic
    [ys, ts] = trap(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,3);
    hold on
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');  
    hold off
    xlabel("Time (t)");
    ylabel("log10(x, y, z)")
    title('ITM of ' + test_name)
    
    % run RK4
    tic
    [ys, ts] = RK4(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,4);
    hold on
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');   
    hold off
    xlabel("Time (t)");
    ylabel("log10(x, y, z)")
    title('RK4 of ' + test_name)
    
    % run AB4 and AM4
    tic
    [ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f);
    toc
    subplot(3,2,5);
    hold on
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');
    hold off
    xlabel("Time (t)");
    ylabel("log10(x, y, z)")
    title('ABM4 of ' + test_name)

end