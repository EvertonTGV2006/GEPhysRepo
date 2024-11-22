%CO02 Gareth Everton

%Ex1: Non linear Equation solver: Lines 3 - 59

tic %start matlab timer

%Mathematical Functions: 
function y = fx(x) %f(x) in Ex1
    y = 3 * exp(-x) -x + 3; % evaluates the function f(x) = 3e^-x - x +3
end

function y = dfx(x) %f'(x) in Ex1
    y = -3 * exp(-x) -1; % evaluates the function f'(x) = df/dx = -3e^-x -1
end

function [xRoot, i] = bisectionRoots(x1, x2, acc) %x1, x2 are inputs of opposite sign, acc is number of dp of accuracy

    %setup inital yValues;
    y1 = fx(x1);
    y2 = fx(x2); 
    i = 0;

    if y2 / y1 > 0
        error("Inputs must evaluate at f(x) to be of opposing sign")
    end
    
    accCheck = false;
    while accCheck ==false
      
        x3 = (x1 + x2) / 2; %evaluate midpoint
    
        y3 = fx(x3); % evaluate new value of f


        if ((y1 > 0 && y3 > 0) || (y1 <0 && y3 <0))
            y1 = y3; % if y1 & y3 are of the same sign, assign y3,x3 to y1,x1
            x1 = x3; % could only reassign x, but copying y eliminates an additional call to f(x) each iteration and is more performant

        else
            y2 = y3; %#ok<NASGU> %if not then must be of same sign as y2 & x2, so assign to those instead
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

        y1 = fx(x1); % calculate values for f(x) and f'(x)
        dy1 = dfx(x1);

        x2 = x1 - (y1 / dy1); % apply newtRaphs method

        i = i+1; %increment iteration counter


        if round(x1, acc) == round(x2, acc)% same acc check as before
            accCheck = true;
        end

        x1 = x2; % set x1 value as input for next iteration
    end

    xRoot = round(x1, acc); % set return value

end


%start of Main Block

fileID = fopen('GE_CO02_Output_EX1.csv', 'w'); %open File for data output


% Find root about x = 3; so take x1 = 2  and x2 = 4, let acc = 10d.p
acc = 10;
[bisRoot, bisIterations] = bisectionRoots(2, 4, acc);
% now take x1 = 3 for newtRaph , acc still = 10;
[newRoot, newIterations] = newtRaph(3, acc);
matRoot = fzero(@fx, 3);

%now write output to file
fprintf(fileID, '%s', ["Method Used:, Calculated Value, Number of Iterations", newline]);
fprintf(fileID, '%s', ["Bisection Method:, ", num2str(bisRoot, 10), ', ', num2str(bisIterations), newline]);
fprintf(fileID, '%s', ["Newton-Raphson Method:, ", num2str(newRoot, 10), ', ', num2str(newIterations), newline]);
fprintf(fileID, '%s', ["MATLAB fzero() Method:, ", num2str(matRoot, 10), newline]);

fclose(fileID);

timeVal = toc;
disp(["Time Elapsed: ", timeVal]) %Print Time Elapsed