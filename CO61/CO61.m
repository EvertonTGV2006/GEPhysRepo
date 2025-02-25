function [ tout , pos, rMin ] = simulate_rocket ( init_pos , init_vel , moon_pos , t)
    % Author: ??? , Date: ??/??/????
    % Simulate the rocket trajectory with the Earth and Moon influence. The coordinate
    % used in this function is centred at Earth’s centre (i.e. Earth centre at (0,0) )
    % and scaled in Moon−radius.
    % The simulation finishes when it simulates for the whole t, or the rocket landed
    % on the Moon.
    % Input:
    % ∗ init_pos: 2−elements vector (x, y) indicating the initial position of the rocket.
    % ∗ init_vel: 2−elements vector (vx, vy) of the initial velocity of the rocket.
    % ∗ moon_pos: a function that receives time, t, and return a 2−elements vector (x, y)
    % (see hint) indicating the Moon position relative to Earth.
    % ∗ t: an N−elements vector of the time step where the position of the rocket will be
    % returned.
    %
    % Output:
    % ∗ tout: an M−elements vector of the time step where the position is described,
    % if the rocket does not land on the Moon, M = N.
    % ∗ pos: (M x 2) matrix indicating the positions of the rocket as function of time,
    % with the first column is x and the second column is y.
    %
    % Example use:
    % >> init_pos = [0, 3.7];
    % >> init_vel = 0.0066 ∗ [cosd(89.9), sind(89.9)];
    % >> moon_pos = @(t) [0, 222];
    % >> t = linspace(0, 10000, 1000);
    % >> [tout, pos] = simulate_rocket(init_pos, init_vel, moon_pos, t);
    % >> plot(

    timesteps = length(t);

    pos(1,:) = init_pos;

    v = init_vel;

    G = 9.63e-7;
    rEarth = [0, 0];
    MEarth = 83.3;
    MMoon = 1;
    
    rMin = 222;

    for i = 2:timesteps
        %compute earth - Sep
        r = pos(i-1,:) - rEarth;
        a = -G * MEarth * r / (sum(r.^2).^(3/2));
        rMoon = moon_pos(t(i));
        r = pos(i-1,:) - rMoon;
        a = -G * MMoon * r / (sum(r.^2).^(3/2)) + a;
        dt = t(i) - t(i-1);
        v = v + a * dt;
        pos(i,:) = pos(i-1,:) + v * dt;
        tout(i) = t(i);
        rMin = min(rMin, sqrt(sum(r.^2)));
        
    end
end

function plotPath(t, p, moon_pos)
    hold off;
    x = p(:,1);
    y = p(:,2);
    col = t;
    scatter(x, y, [], col, 'fill');
    hold on;



    %plot moon;
    mp = zeros(length(t),2);

    for i = 1:length(t)
        mp(i,:) = moon_pos(t(i));
    end

    mx = mp(:,1);
    my = mp(:,2);
    scatter(mx, my, [], col, 'fill');

    axis([-300 300 -300 300])
    shg;
    grid on;
end

moon_pos = @(t) [222*cos(2.6615e-6*t), 222*sin(2.6615e-6*t)];
% moon_pos = @(t) [0, 222];

function [t, p, theta] = findLaunchAngle(moon_pos)
    v0 = 0.0066;
    thetas = 0:pi/4:pi;
    rMin = 222;
    thetaRes = 0;
    for i = 1:length(thetas)
        v = [v0 * cos(thetas(i)), v0 * sin(thetas(i))];
        [t, p, rMr] = simulate_rocket([0, 3.7], v, moon_pos, 0:10:2000000);
        plotPath(t, p, moon_pos);
        thetas(i)
        rMr
        if rMr < rMin
            thetaRes = thetas(i);
            rMin = rMr;
        end
    end
    dtheta = 0.1
    theta2 = thetaRes + dtheta;
    drdt = 0;
    for i = 1:10
        v = [v0 * cos(theta2), v0 * sin(theta2)];
        [t, p, rMr2] = simulate_rocket([0, 3.7], v, moon_pos, 0:10:2000000);
        plotPath(t, p, moon_pos);
        drdt = (rMr2 - rMr )/dtheta
        dtheta = -rMr2 / drdt
        if (dtheta > 0.1)
            dtheta = 0.1
            dtheta = dtheta * rMr2;
        end
        if (dtheta < -0.1)
            dtheta = -0.1
            dtheta = dtheta * rMr2;
        end
        
        theta2 = theta2 + dtheta
        rMr = rMr2
        if(rMr < 0.1)
            break
        end
    end
    theta = theta2;

end


[t, p, theta] = findLaunchAngle(moon_pos);

"Best launch angle: "
theta_deg = theta * 180 / pi


