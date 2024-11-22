%CO02 Gareth Everton

tic %Start MATLAB timer

%Ex2: Simpsons Rule


function y = gx(x) % erf(x) for Ex2
    y = exp(-(x.^2)) * 2 * pi.^-0.5;
end
    
function y = hx(x) % 2nd func for Ex2
    y = (sin(100 * pi * x)).^2 * exp(-(x.^2));
    y = y * 2 * pi.^-0.5;
end

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

fileID = fopen('GE_CO02_Output_EX2.csv', 'w'); %open File for data output
    
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

fclose(fileID);

timeVal = toc;
disp(["Time Elapsed: ", timeVal]) %Print Time Elapsed
