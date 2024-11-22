%CO02 Gareth Everton

tic %Start MATLAB timer

%Exercise 3: Curve Fitting


function y = lin(x, a, b)
    y = a + b * x;
end

function [a, b] = linearCurveFit(x, y)
    N = length(x);
    if length(y) ~= N
        error("Input arrays of unequal length");
    end

    %zero - intialize variables
    sumX = 0;
    sumY = 0;
    sumXsq = 0;
    sumXY = 0;

    %calculate sums over the ranges
    for i = 1: N
        sumX =sumX + x(i);
        sumY = sumY + y(i);
        sumXsq = sumXsq + x(i).^2;
        sumXY = sumXY + (x(i) * y(i));
    end

    %create the matrix from the required sums
    matX = [N, sumX; sumX, sumXsq];
    vecY = [sumY; sumXY];

    vecA = matX\vecY; %Perform matrix back multiplation using matlab \ operator

    %assign return values from answer vector
    a = vecA(1);
    b = vecA(2);
end

%helper function to calcualte rms errors
function r = calculateRMS(y, yFit)

    rSq = 0; %zero intialize variables

    for i = 1 : length(y)
    rSq = rSq + (yFit(i) - y(i)).^2; %perform the sum squares
    end

    r = rSq.^0.5; %root the return value at the end
    r = r /length(y); %divide by number of points N
end


function linearFitExecution(fileID, x, y)

    %run linear fit on x and y values
    [a, b] = linearCurveFit(x, y);

    yFit = zeros(1, length(x)); %Pre allocate for memory performance

    for i = 1 : length(x)
        yFit(i) = lin(x(i), a, b); %Use helper function lin to create fitted values for the straight line
    end

    r = calculateRMS(y, yFit); %calculate final RMS

    figure();
    plot(x, y, '.b'); %now plot data to visualise
    hold on;
    plot(x, yFit, 'r-')
    annotation('textbox',[.2 .9 0 0],'String',["RMS for Line:", r],'FitBoxToText','on');
    annotation('textbox',[.6 .4 0 0],'String',["y = ",a," + ",b,"x"], 'FitBoxToText','on');

    fprintf(fileID, '%s', ['X, Y, Y_FIT, r = ',num2str(r),', a = ',num2str(a), ', b = ',num2str(b), ',', newline]); %and write out raw Data to .csv
    for i = 1 : length(x)
        fprintf(fileID, '%s', [num2str(x(i)), ', ', num2str(y(i)), ', ', num2str(yFit(i)), newline]);
    end
end


function [x, y] = polynomialGenerator(N, m, a, e)
    if m + 1 ~= length(a)
        %error(["Incorrect number of coefficients: Polynomial of order ",num2str(m)," was specified but ",num2str(length(a))," coefficients were specified"]);
        error("Incorrect number of coefficients");
    end

    x = zeros(1, N+1); %initialize arrays for memory performance
    y = zeros(1, N+1);

    for i = 0 : N
        x(i + 1) = -5  + 10 * i / N; %Calculate our x value using linearInterpolation
        y(i + 1) = 0; %initialise y value to 0
        for j  = 1 : length(a)
            y(i + 1) = y(i +1) + a(j) * x(i + 1).^(j-1); %Sum over coefficients adding a * x ^ (j-1) to the y value
        end
        %this gives us an exact fitted polynomial

        y(i + 1) = y(i + 1) + e * randn(); % Adds noise to the final value
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
        for j = 1 : length(x) %calculate all the sums we need for the matrix MatX
            ySums(i) = ySums(i) + y(j) * x(j).^(i-1);
            xSums(i) = xSums(i) + x(j).^(i-1); %Perform xSum within one loop to reduce loop overhead and make use of AVX (hopefully - not sure if it makes a difference in MATLAB but defo does in C++)
            xSums(i + n) = xSums(i + n) + x(j).^(i-1+n);
        end
    end

    for i = 1:n
        for j = 1:n
            matX(i, j) = xSums(i+j-1); %Assign values to the elements of matX;
        end
    end

    a = matX\ySums';%Perform back multiplication on the transpose of ySums to find a, our vector of coefficients

    for i = 1 : length(x)
        for j  = 1 : n
            yFit(i) = yFit(i) + a(j) * x(i).^(j-1); %now using our vector of coefficients we calculate fitted values of y to feed into our plot
        end
    end

    figure()
    plot(x, y, 'b.'); %plot to visulaise data
    hold on
    plot(x, yFit, 'r-');
    r = calculateRMS(y, yFit);
    annotation('textbox',[.2 .9 0 0],'String',["RMS for Polynomial:", r],'FitBoxToText','on')
    annotation('textbox',[.6 .9 0 0],'String',["y = ",a(1)," + ",a(2),"x + ",a(3),"x^2 + ",a(4),"x^3"], 'FitBoxToText','on');
    


end

function polynomialFitExecute(fileID)
    %function to call the above polynomial fit functions
    aReal = [1.5, -2.5, 0.7, -1.2]; %Initialise our actual aReal coefficients

    [x, y] = polynomialGenerator(100, 3, aReal, 1); % generate a polynomial with noise based on these coefficients

    [aFit, yFit, r] = polynomialFit(x, y, 3); %Run the polynomial fitting function

    %now write Data out to the .csv
    fprintf(fileID, '%s', 'X:,Y:,Y_FIT,y = ');
    for i = 1 : length(aFit)
        fprintf(fileID, '%s', [num2str(aFit(i), 3), +'x^', num2str(i-1, 1), ' ']);
    end
    fprintf(fileID, '%s', [',r = ',num2str(r, 5), newline]);
    for i = 1 : length(x)
        fprintf(fileID, '%s', [num2str(x(i), 3),', ',num2str(y(i), 3),', ',num2str(yFit(i), 3),newline]);
    end

end

%start of Main block

fileID = fopen('GE_CO02_Output_EX3.csv', 'w'); %open File for data output

data = load("linear.csv"); %Load sample input data for the linear case

x = data(:, 1); %assign the linear data to x and y variables
y = data(:, 2);

linearFitExecution(fileID, x, y); %run the function to fit the curve

%Now do polynomial fitting

polynomialFitExecute(fileID); % run the function to fit the curve

fclose(fileID);

timeVal = toc;
disp(["Time Elapsed: ", timeVal]) %Print Time Elapsed




