axis off;

% define pixel dimensions
width  = 40;
height = 40;

% timestep and time variables
h = 1e-6;
t_0 = 0;
t_max = 8;

% 
f = @o;
y_0 = [0 0.23 0.0]';

petri_dish = zeros(width, height, 3);
petri_dish(floor(width / 2), floor(height/2), :) = y_0;

diffusion_kernel = [
    1./80. 1./80. 1./80.;
    1./80. 9./10. 1./80.;
    1./80. 1./80. 1./80.
];

diffusion_kernel2 = [
    0.07 0.1  0.07
    0.1  0.32 0.1
    0.07 0.1  0.07
];

ts = (t_0:h:t_max)';

imwrite(petri_dish, ['oregonator_step_0', '.png']);

for k = 1:length(ts)-1
    % make a frame every 100000 steps
    if mod(k, 100000) == 0
        disp(k);
        petri_dish(:, :, 1) = conv2(petri_dish(:, :, 1), diffusion_kernel2, 'same');
        petri_dish(:, :, 2) = conv2(petri_dish(:, :, 2), diffusion_kernel2, 'same');
        petri_dish(:, :, 3) = conv2(petri_dish(:, :, 3), diffusion_kernel2, 'same');
        imwrite(petri_dish, [sprintf('oregonator_step_%d', k), '.png']);
    end

    for i = 1:width
        for j = 1:height
            % Perform Euler Step
            next = h*f(ts(k), [petri_dish(i, j, 1) petri_dish(i, j, 2) petri_dish(i, j, 3)]);

            petri_dish(i, j, 1) = petri_dish(i, j, 1) + next(1);
            petri_dish(i, j, 2) = petri_dish(i, j, 2) + next(2);
            petri_dish(i, j, 3) = petri_dish(i, j, 3) + next(3);
    
            % Perform convolution for each concentration
        end
    end
end

imwrite(petri_dish, [sprintf('oregonator_step_%d', k), '.png']);

% plot_test(t_0, y_0, h, t_max, f, "BZ Model");


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
function [] = plot_test(t_0, y_0, h, t_max, f, test_name)
    clf;
    hold on;

    % run FEM
    [ys, ts] = FEM(t_0, y_0, h, t_max, f);
    % subplot(3,2,1);
    plot(ts, ys(1, :), 'r.-');
    % plot(ts, ys(2, :), 'g.-');
    plot(ts, ys(3, :), 'b.-');
    title('FEM of ' + test_name)
end