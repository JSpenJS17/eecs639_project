h = .000005;
t_0 = 0;
t_max = 25;
f = @dyz;
y_0 = [1; 1; 1;];
plot_test(t_0, y_0, h, t_max, f, "BZ Model");


% plot concentration of Y with respect to X
% plot(ys(:,1), ys(:,2), 'green');

function [dydt] = dyz(t, y)
    %this function works with these starting parameters, but it's using the
    %scaled version instead of the regular one given in the project

    % reaction parameters
    a = 0.06;
    b = 0.02;

    kc = 1.0;
    k2 = 2.4E+06;
    k3 = 1.28;
    k4 = 3.0E+03;
    k5 = 33.6;

    %calculated scaling stuff
    eta1 = kc * b / k5 / a;
    eta2 = 2.0 * kc * k4 * b / k2 / k5 / a;
    q = 2.0 * k3 * k4 / k2 / k5;
    f = 1.0;

    % grab x, y, and z
    u = y(1);
    v = y(2);
    w = y(3);
    
    % calculate new x, y, and z
    dudt = (   q * v - u * v + u * ( 1.0 - u ) ) / eta1;
    dvdt = ( - q * v - u * v + f * w ) / eta2;
    dwdt = u - w;
    
    % return them
    dydt = [ dudt; dvdt; dwdt ];

end

% define helper function to run solvers and plot results
function [] = plot_test(t_0, y_0, h, t_max, f, test_name)
    figure;

    % run FEM
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    subplot(3,2,1);
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');
    title('FEM of ' + test_name)
    
    % run BEM,  disabled until root finding is finished
    [ys, ts] = BEM(t_0, y_0, h, t_max, f);
    subplot(3,2,2);
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');
    title('BEM of ' + test_name)
  
    % run trap on electric circuit, disabled until fixed
    [ys, ts] = trap(t_0, y_0, h, t_max, f);
    subplot(3,2,3);
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');  
    title('Trap of ' + test_name)
    
    % run RK4
    [ys, ts] = RK4(t_0, y_0, h, t_max, f);
    subplot(3,2,4);
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');   
    title('RK4 of ' + test_name)
    
    % run AB4 and AM4
    [ys, ts] = predictor_corrector(t_0, y_0, h, t_max, f);
    subplot(3,2,5);
    plot(ts, log10(ys(1, :)), 'r-');
    plot(ts, log10(ys(2, :)), 'g-');
    plot(ts, log10(ys(3, :)), 'b-');  
    title('PC of ' + test_name)

end