x = [300.4, 280.2, 260.2, 240.0, 220.1, 200];
xerr = 0.05 * ones(length(x));
y = [3.413, 3.341, 3.13, 2.832, 2.651, 2.424];
yerr = 0.0005 * ones(length(y));

I_2 = y.^2;
I_2err = 2 * (yerr / y) * I_2;

errorbar(x, I_2, I_2err, I_2err, xerr, xerr);

c = polyfit(x, I_2, 1);
xVal = [0, 100, 200, 300];
yFit = c(1) * xVal + c(2);

hold on;

plot(xVal, yFit);