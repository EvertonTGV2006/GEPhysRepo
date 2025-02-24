function [ tout , pos ] = simulate_rocket ( init_pos , init_vel , moon_pos , t)
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

    end


end


moon_pos = @(t) [222*cos(2.6615e-6*t), 222*sin(2.6615e-6*t)];

v0 = 0.0066;
theta = 89.9 * pi/180;
v = [v0 * cos(theta), v0 * sin(theta)];
[t, p] = simulate_rocket([0, 3.7], v, moon_pos, 0:10:200000);
plot(p(:,1),p(:,2))
axis([-300 300 -300 300])



    
