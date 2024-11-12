%CO02 Gareth Everton

%Ex1: Non linear Equation solver: Lines 3 - 59
%Ex2: Numerical Integration: Lines 60 - 




%Mathematical Functions: 
function y = fx(x) %f(x) in Ex1
    y = 3 * exp(-x) -x + 3; % evaluates the function f(x) = 3e^-x - x +3
end

function y = dfx(x) %f'(x) in Ex1
    y = -3 * exp(-x) -1; % evaluates the function f'(x) = df/dx = -3e^-x -1
end

function y = gx(x) % erf(x) for Ex2
y = exp(-(x.^2)) * 2 * pi.^-0.5;
end

function y = hx(x) % 2nd func for Ex2
y = (sin(100 * pi * x)).^2 * exp(-(x.^2));
y = y * 2 * pi.^-0.5;
end

function y = lin(x, a, b)
    y = a + b * x;
end




function [xRoot, i] = bisectionRoots(x1, x2, acc) %x1, x2 are inputs of opposite sign, acc is number of dp of accuracy
y1 = fx(x1);
y2 = fx(x2);
i = 0;

accCheck = false;
while accCheck ==false
  
    x3 = (x1 + x2) / 2; %evaluate midpoint

    y3 = fx(x3); % evaluate new value of f
    if ((y1 > 0 && y3 > 0) || (y1 <0 && y3 <0))
        y1 = y3; % if y1 & y3 are of the same sign, assign y3,x3 to y1,x1
        x1 = x3; % could only reassign x, but copying y eliminates an additional call to f(x) each iteration and is more performant
    else
        y2 = y3; %if not then must be of same sign as y2 & x2, so assign to those instead
        x2 = x3;
    end
    if round(x1, acc) == round(x2, acc)% now carry out accuracy checks, if the rounded versions match we have converged on a root to acc many d.p. so we can break out of the loop
        accCheck = true;
    end
    i = i+1;
end
xRoot = round(x1, acc); % set return value
end

function [xRoot, i] = newtRaph(x1, acc) %similar to before, but 1 less input

i = 0;
accCheck = false;
while accCheck == false
    y1 = fx(x1); % calculate values
    dy1 = dfx(x1);
    x2 = x1 - (y1 / dy1); % apply newtRaphs method
    i = i+1; %increment iteration counter
    if round(x1, acc) == round(x2, acc)% same acc check as before
        accCheck = true;
    end
    x1 = x2;
end
xRoot = round(x1, acc); % set return value
end



%Ex2: Simpsons Rule

function area = simpsonsRule(fx, x1, x2, n)
h = (x2 - x1) / n;
y1 = fx(x1); %evaluate first and last values
y2 = fx(x2);
workingSum = y1 + y2; %start working sums %"ends x 1"
for i = 1:(n-1)
    if mod(i, 2) == 0 
        workingSum = workingSum + 2 * fx(x1 + i * h); %"evens" * 2
    end
    if mod(i,2) == 1
        workingSum = workingSum + 4 * fx(x1 + i * h); %"odds"  * 4
    end
end
area = workingSum * h / 3;
end




%Ex 3

function [a, b] = linearCurveFit(x, y)
    N = length(x);
    if length(y) ~= N
        error("Input arrays of unequal length");
    end

    sumX = 0;
    sumY = 0;
    sumXsq = 0;
    sumXY = 0;

    for i = 1: N
        sumX =sumX + x(i);
        sumY = sumY + y(i);
        sumXsq = sumXsq + x(i).^2;
        sumXY = sumXY + (x(i) * y(i));
    end

    matX = [N, sumX; sumX, sumXsq];
    vecY = [sumY; sumXY];

    vecA = matX\vecY;
    a = vecA(1);
    b = vecA(2);
end

function r = calculateRMS(y, yFit)
    rSq = 0;
    for i = 1 : length(y)
    rSq = rSq + (yFit(i) - y(i)).^2;
    end
    r = rSq.^0.5;
end


function linearFitExecution(fileID, x, y)
    [a, b] = linearCurveFit(x, y);
    yFit = zeros(1, length(x)); %Pre allocate for memory performance
    for i = 1 : length(x)
        yFit(i) = lin(x(i), a, b);
    end
    r = calculateRMS(y, yFit);

    plot(x, y, '.b');
    hold on;
    plot(x, yFit, 'r-')

    fprintf(fileID, '%s', ['X, Y, Y_FIT, r = ',num2str(r),', a = ',num2str(a), ', b = ',num2str(b), ',', newline]);
    for i = 1 : length(x)
        fprintf(fileID, '%s', [num2str(x(i)), ', ', num2str(y(i)), ', ', num2str(yFit(i)), newline]);
    end
end


function [x, y] = polynomialGenerator(N, m, a, e)
    if m + 1 ~= length(a)
        %error(["Incorrect number of coefficients: Polynomial of order ",num2str(m)," was specified but ",num2str(length(a))," coefficients were specified"]);
        error("Incorrect number of coefficients");
    end
    x = zeros(1, N+1);
    y = zeros(1, N+1);
    for i = 0 : N
        x(i + 1) = -5  + 10 * i / N;
        y(i + 1) = 0;
        for j  = 1 : length(a)
            y(i + 1) = y(i +1) + a(j) * x(i + 1).^(j-1);
        end
        y(i + 1) = y(i + 1) + e * randn();
    end
end

function [a, yFit, r] = polynomialFit(x, y, m)
    %Want to minimise compute time for these matrices, so preallocate arrays
    n = m + 1; %Polynomial order m has n cofficients


    xSums = zeros(1, 2 * n);
    ySums = zeros(1, n);
    yFit = zeros(1,length(x));
    matX = zeros(n, n);
    for i = 1 : n
        for j = 1 : length(x)
            ySums(i) = ySums(i) + y(j) * x(j).^(i-1);
            xSums(i) = xSums(i) + x(j).^(i-1); %Perform xSum within one loop to reduce loop overhead and make use of AVX (hopefully - not sure if it makes a difference in MATLAB but defo does in C++)
            xSums(i + n) = xSums(i + n) + x(j).^(i-1+n);
        end
    end
    for i = 1:n
        for j = 1:n
            matX(i, j) = xSums(i+j-1);
        end
    end

    a = matX\ySums';

    for i = 1 : length(x)
        for j  = 1 : n
            yFit(i) = yFit(i) + a(j) * x(i).^(j-1);
        end
    end
    plot(x, y, 'b.');
    hold on
    plot(x, yFit, 'r-')

    r = calculateRMS(y, yFit);
    


end

function polynomialFitExecute(fileID)
    aReal = [1.5, -2.5, 0.7, -1.2];
    [x, y] = polynomialGenerator(100, 3, aReal, 1);
    [aFit, yFit, r] = polynomialFit(x, y, 3);
    fprintf(fileID, '%s', 'X:,Y:,Y_FIT,y = ');
    for i = 1 : length(aFit)
        fprintf(fileID, '%s', [num2str(aFit(i), 3), +'x^', num2str(i-1, 1), ' ']);
    end
    fprintf(fileID, '%s', [',r = ',num2str(r, 5), newline]);
    for i = 1 : length(x)
        fprintf(fileID, '%s', [num2str(x(i), 3),', ',num2str(y(i), 3),', ',num2str(yFit(i), 3),newline]);
    end

end
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










fileID = fopen('GE_CO02_Output.csv', 'w'); %open File for data output
%EXERCISE_NUMBER = input("Select Exercise: (1 - 4): ");
EXERCISE_NUMBER = 4;
%Main code
switch EXERCISE_NUMBER
    case 1
        % Find root about x = 3; so take x1 = 2  and x2 = 4, let acc = 3d.p
        acc = 10;
        [bisRoot, bisIterations] = bisectionRoots(2, 4, acc);
        % now take x1 = 3 for newtRaph , acc still = 3;
        [newRoot, newIterations] = newtRaph(3, acc);
        matRoot = fzero(@fx, 3);
    
        %now write output to file
        fprintf(fileID, '%s', ["Method Used:, Calculated Value, Number of Iterations", newline]);
        fprintf(fileID, '%s', ["Bisection Method:, ", num2str(bisRoot, 10), ', ', num2str(bisIterations), newline]);
        fprintf(fileID, '%s', ["Newton-Raphson Method:, ", num2str(newRoot, 10), ', ', num2str(newIterations), newline]);
        fprintf(fileID, '%s', ["MATLAB fzero() Method:, ", num2str(matRoot, 10), newline]);


    case 2     %Ex 2, simpsons rule

    %Initialse all the variables
    X_Values = transpose([0.2, 0.4, 0.6, 0.8, 1.0]);
    N_1_Values = zeros(length(X_Values),1);
    N_10_Values = zeros(length(X_Values),1);
    N_100_Values = zeros(length(X_Values),1);
    ERF_Values = zeros(length(X_Values),1);

    for iX = 1:length(X_Values) %Calculate area for each x input using simpsons rule and store it in an array
        N_1_Values(iX,1) = simpsonsRule(@gx, 0, X_Values(iX), 1);
        N_10_Values(iX,1) = simpsonsRule(@gx, 0, X_Values(iX), 10);
        N_100_Values(iX,1) = simpsonsRule(@gx, 0, X_Values(iX), 100);
        ERF_Values(iX,1) = erf(X_Values(iX));
    end
    %write output to csv
    fprintf(fileID, '%s',["Integral of Error Function (Part A) using Simpson's Rule: X, N = 1, N = 10, N = 100, MATLAB erf(x)", newline]);
    for iX = 1:length(X_Values)
        fprintf(fileID, '%s', [num2str(X_Values(iX), 2), ', ', num2str(N_1_Values(iX), 10), ', ', num2str(N_10_Values(iX), 10), ', ', num2str(N_100_Values(iX), 10), ', ', num2str(ERF_Values(iX), 10), newline]);
    end
    fprintf(fileID, '%s', newline);

    % now calculate for oscillatory function
    for iX = 1:length(X_Values)
        N_1_Values(iX,1) = simpsonsRule(@hx, 0, X_Values(iX), 1);
        N_10_Values(iX,1) = simpsonsRule(@hx, 0, X_Values(iX), 10);
        N_100_Values(iX,1) = simpsonsRule(@hx, 0, X_Values(iX), 100);
    end
    %write output to csv
    fprintf(fileID, '%s', ["Integral of Oscillatory Function (Part B) using Simpon's rule: X, N = 1, N = 10, N = 100", newline]);
    for iX = 1:length(X_Values)
        fprintf(fileID, '%s', [num2str(X_Values(iX), 2), ', ', num2str(N_1_Values(iX), 10), ', ', num2str(N_10_Values(iX), 10), ', ', num2str(N_100_Values(iX), 10), newline]);
    end

    case 3
        %Exercise 3: Curve Fitting
        %CALC_ORDER = input("Order of the best fit line to be calculated: 1 = Linear, 2 = Quadratic etc: ");
        CALC_ORDER = 2;
        switch CALC_ORDER
            case 1
                data = load("linear.csv");

                x = data(:, 1);
                y = data(:, 2);

                linearFitExecution(fileID, x, y);
            case 2
                polynomialFitExecute(fileID);


        end

    case 4
        %HALF_LIFE_X = input("Half life of element X (Years): ");
        %HALF_LIFE_Y = input("Half life of element Y (Years): ");

        HALF_LIFE_X = 2;
        HALF_LIFE_Y = 16;

        [t, x, y, z] = executeNuclearDecay(50, 1, HALF_LIFE_X, HALF_LIFE_Y, 1000);
        fprintf(fileID, '%s', ["t, x, y, z",newline]);
        for i = 1 : length(t)
            fprintf(fileID, '%s', [num2str(t(i)), ', ', num2str(x(i), 3), ', ', num2str(y(i), 3), ', ', num2str(z(i), 3), newline]);
        end

        plot(t, x, 'b-');
        hold on;
        plot(t, y, 'r-');
        plot(t, z, 'g-');


end
fclose(fileID);



