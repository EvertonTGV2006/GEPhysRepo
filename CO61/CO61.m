function [ tout , pos, rMin ] = simulate_rocket ( init_pos , init_vel , moon_pos , t)
    % Author: Gareth , Date: 28/02/25
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

    tic %Start iteration timer

    timesteps = length(t); %Number of timesteps is the count of the number of times given in array t;

    pos(1,:) = init_pos; %set starting position in p array

    v_0 = init_vel; %Set v_0 as inital velocity, initialise values for euler integrations

    G = 9.63e-7; %Constant value of G
    rEarth = [0, 0]; %Fix Position of Earth to be 0, ingnore effect of moon's gravity on Earth
    MEarth = 83.3; %Mass of the Earth set to be 83.3
    MMoon = 1; %Mass of the moon set to be 1;
    
    rMin = 222; %Set large inital value of rMin so that it can be beaten by better orbits

    tout = zeros(size(t)); %Initialize array of t values as 0s

    for i = 2:timesteps
        
        p_0 = pos(i-1,:); %Starting position for iteration

        %First Euler iteration
        r_1 = p_0 - rEarth; %Find Earth-Sat separation
        a_1 = -G * MEarth * r_1 / (sum(r_1.^2).^(3/2)); %Compute acceleration due to Earth
        rMoon = moon_pos(t(i-1)); %Get Moon Position in space for given time
        r_1 = p_0 - rMoon; %Get Moon-Sat Separation
        a_1 = -G * MMoon * r_1 / (sum(r_1.^2).^(3/2)) + a_1; %Add acceleration due to Moon
        dt = t(i) - t(i-1); %Compute timestep dt
        dv_1 = a_1 * dt; %find dv_1
        v_1 = v_0 + dv_1; %Increment velocity
        dp_1 = v_1 * dt; %Find dp
        p_1 = p_0 + dp_1; %Increment Positions

        %2nd Euler iteration, same code as above but using v_1 and p_1 and t(i) instead
        r_2 = p_1 - rEarth;
        a_2 = -G * MEarth * r_2 / (sum(r_2.^2).^(3/2));
        rMoon = moon_pos(t(i));
        r_2 = p_1 - rMoon;
        a_2 = -G * MMoon * r_2 / (sum(r_2.^2).^(3/2)) + a_2;
        dt = t(i) - t(i-1);
        dv_2 = a_2 * dt;
        v_2 = v_0 + dv_2;
        dp_2 = v_2 * dt;


        dv_t = 0.5* (dv_1 + dv_2); %Combine dv_x into total dv_t
        dp_t = 0.5* (dp_1 + dp_2); %Same for position

        p_t = p_0 + dp_t; %Update position
        v_t = v_0 + dv_t; %Update velocity


        v_0 = v_t; %Write new velocity to velocity for next iteration
        pos(i,:) = p_t; %Store position

        tout(i) = t(i); %Write out timestep
        rMin = min(rMin, sqrt(sum(r_2.^2))); %Compute rMin for return value
        
    end

    toc %Output elapsed time for iteration
end

function plotPath(t, p, moon_pos)
     %Author: Gareth     Date: 28/02/25
     %Input Paramaters:
     %t: M size array of double values indicating t values at positions descriped in p
     %p: M x 2 size array of double values indicating x and y positions in space at time t
     %moon_pos: function handle to function that takes input value of t and returns x and y coordinates of the moon.
    
    hold off; %Hold off to clear screen for new plot
    x = p(:,1); %Initialize x and y arrays for plot
    y = p(:,2);
    col = t; %Set colour of line to time at that position
    scatter(x, y, [], col, 'fill'); %plot scatter of positions with colour gradient showing dt for input position
    hold on;
    %Keep plot on screen


    %plot moon;
    mp = zeros(length(t),2); %Initialize array for position of moon
    for i = 1:length(t)
        mp(i,:) = moon_pos(t(i)); %Populate moon position vector mp with values for given postions at time t
    end

    mx = mp(:,1); %Initialsie x and y values for the moon, col still the same as used above
    my = mp(:,2);
    scatter(mx, my, [], col, 'fill'); %Same scatter plot but for moon

    axis([-300 300 -300 300]) %Set axis of plot so that it is scaled properly
    c = colorbar; %Add a colorbar to show timescale of orbits
    c.Label.String = 'Time / s'; %Label colorbar with Time values
    grid on; %Set the grid to be on
    shg; %Show the completed figure to the user.
end

moon_pos = @(t) [222*cos(2.6615e-6*t), 222*sin(2.6615e-6*t)];
% moon_pos = @(t) [0, 222];

function [t, p, theta] = findLaunchAngle(moon_pos)
    %Author: Gareth     Date: 28/02/25
    %Input parameters: moon_pos: handle to function that returns a 2D position vector for the moon for a given input time t;
    v0 = 0.0066; %Initial velocity
    thetas = 0:pi/4:pi; %Seed initial "guesses" for theta
    rMin = 222; %Set rMin closest approach to start value
    thetaRes = 0; %resultant theta
    for i = 1:length(thetas)
        v = [v0 * cos(thetas(i)), v0 * sin(thetas(i))]; %Convert v0 and theta to vector
        [t, p, rMr] = simulate_rocket([0, 3.7], v, moon_pos, 0:10:2000000); %Simulate trajectory using simulate_rocket function as above
        plotPath(t, p, moon_pos); %Plot the resultant path so that the user can visualise the currently calculated trajectory;
        thetas(i) %Print out current theta
        rMr %Print out minimun radius obtained from the current simulation iteration
        if rMr < rMin %If this closes approach is better than previous ones, save the value of theta for use later, and update new closes approach
            thetaRes = thetas(i); 
            rMin = rMr;
        end
    end
    dtheta = 0.1 %Print seed value of dtheta
    theta2 = thetaRes + dtheta; %Now next guess of theta is based on previous guess and dtheta
    drdt = 0; %Rate of change of r wrt theta, where r is closest approach to moon.
    for i = 1:10 % set a limit of 10 iterations
        v = [v0 * cos(theta2), v0 * sin(theta2)]; %Calculate seed velocity from v0 and theta
        [t, p, rMr2] = simulate_rocket([0, 3.7], v, moon_pos, 0:10:2000000); %simulate trajectory
        plotPath(t, p, moon_pos); %Plot path for user
        drdt = (rMr2 - rMr )/dtheta %Print rate of change of closest approac wrt theta for user information
        dtheta = -rMr2 / drdt %calculate new dtheta for drdt
        if (dtheta > 0.1)
            dtheta = 0.1
            dtheta = dtheta * rMr2; %Do some bounds checks so that dtheta is not too large, this adds stability to the simulation
        end
        if (dtheta < -0.1)
            dtheta = -0.1
            dtheta = dtheta * rMr2;
        end
        
        theta2 = theta2 + dtheta %Print new theta
        rMr = rMr2 %Update disrtance of closes approach
        if(rMr < 0.1) %Check if we have met our success criteria of being within 0.1 radii of the moon.
            break
        end
    end
    theta = theta2; %Return best value of theta that was found

end


[t, p, theta] = findLaunchAngle(moon_pos); %Find best launch angle for rocket.

"Best launch angle: " %Print out best launch angle
theta_deg = theta * 180 / pi


