t = readtable("3_1.txt", "FileType","text",'Delimiter', '\t')

peak = table2array(t(:,"Peak"));
time = table2array(t(:, "Time"));
angle = table2array(t(:,"Angle"));


peak_even = peak(2:2:end);
peak_odd = peak(1:2:end);
time_even = time(2:2:end);
time_odd = time(1:2:end);
angle_even = angle(2:2:end);
angle_odd = angle(1:2:end);

function y = expoff(x,a,b,c)
    y = a * pow(b, x) + c;
end

bestfit = fit(peak_even, angle_even, @expoff)