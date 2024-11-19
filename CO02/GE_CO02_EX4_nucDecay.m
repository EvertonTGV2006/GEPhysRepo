%%CO02 Gareth Everton

tic %Start MATLAB timer

%Exercise 4: Nulcear Decay

function [x, y, z] = evaluateDecayEvent(t, aX, aY, x0)
    x = x0 * exp(-aX * t);
    if aX == aY
        y = x0 * aY * t * exp(-aY * t);
    else
        y = x0 * aX * ((exp(-aX * t) - exp(-aY * t)) / (aY - aX));
    end
    z = x0 - x- y;
end

function [t, x, y, z]  = executeNuclearDecay(tF, dt, tX, tY, x0)
    aX = log(2)/tX;
    aY = log(2)/tY;
    count = ceil(tF/dt) + 1;
    t = zeros(1,count);
    x = zeros(1,count);
    y = zeros(1,count);
    z = zeros(1,count);
    x(1) = x0;
    for i = 1 : count
        [x(i+1), y(i+1), z(i+1)] = evaluateDecayEvent(i * dt, aX, aY, x0);
        t(i+1) = i * dt;
    end
end

function [x, y, z] = nuclearDecayStep(x, y, z, dt, ax, ay)
    px = 1 - exp(-ax * dt);
    py = 1 - exp(-ay * dt);
    dx = 0;
    dy = 0;
    for i = 1 : x
        if px > rand()
            dx = dx + 1;
        end
    end
    for j = 1 : y
        if py > rand()
            dy = dy + 1;
        end
    end
    z = z + dy;
    y = y + dx - dy;
    x = x - dx;
end


%Start of Main block 

fileID = fopen('GE_CO02_Output.csv', 'w'); %open File for data output

%HALF_LIFE_X = input("Half life of element X (Years): ");
%HALF_LIFE_Y = input("Half life of element Y (Years): ");

HALF_LIFE_X = 2;
HALF_LIFE_Y = 16;

[t, x, y, z] = executeNuclearDecay(50, 1, HALF_LIFE_X, HALF_LIFE_Y, 1000); %Run the first nuclear decay with analytical solutions

%Output Data to file
fprintf(fileID, '%s', ["t, x, y, z",newline]);
for i = 1 : length(t)
    fprintf(fileID, '%s', [num2str(t(i)), ', ', num2str(x(i), 3), ', ', num2str(y(i), 3), ', ', num2str(z(i), 3), newline]);
end


%Plot X,Y,Z on a graph to visualise data for the analytical solution
plot(t, x, 'b-');
hold on;
plot(t, y, 'r-');
plot(t, z, 'g-');
hold off;

pause(10);
%Now for the animation / time step probabilistic approach


%Initialise variables
ax = log(2) / HALF_LIFE_X;
ay = log(2) / HALF_LIFE_Y;
STEP = 0.1;
N = 50;
x = 1000;
y = 0;
z = 0;

%And execute loop N / step times, in this case 500 times
for i = 1 : N / STEP
    [x, y, z] = nuclearDecayStep(x, y, z, STEP, ax, ay); %Run the decay probabilites on current population values;

    %plot on a graph
    bar(i * STEP, [x, y, z], "stacked");

    xleg = "X: " + num2str(x);
    yleg = "Y: " + num2str(y);
    zleg = "Z: " + num2str(z);
    legend(xleg, yleg, zleg);

    %Add the pause to control animation speed
    pause(0.05);
end


fclose(fileID);

timeVal = toc;
disp(["Time Elapsed: ", timeVal]) %Print Time Elapsed



