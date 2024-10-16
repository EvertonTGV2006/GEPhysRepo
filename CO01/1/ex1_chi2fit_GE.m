% CO01 - Exercise 1 - Chi-Squared

close all;

%Load Data
myData = load("chi2fit.csv");

%Assign Data to vectors;
vecX = myData(:,1);
vecY = myData(:,2);
vecF = myData(:,3);

% set delta_y
delta_y = 5;

%Initialize variable Chi Squared and error variables;
ChiSquared = 0;
error = delta_y * ones(1,length(vecX)); %assumes size of error bar = delta_y

%iterate over each element, with i being the index, and add to chiSquared sum
for i = 1:length(vecX)
    ChiSquared =ChiSquared + (((vecY(i)-vecF(i)).^2)/(delta_y.^2));
end


fig1 = figure(); %initialize figure
errorbar(vecX, vecY, error, 'vertical', '.'); %plot raw Y as points with errorBars
hold on; %hold figure for next plot 
plot(vecX, vecF, 'r-'); %plot fitted line on graph

% Labelling Graph
xlabel('X Value'); %Add X Label
ylabel('Y Value'); %Add Y Label
xlim([-5 5]);
ylim([-50 20]);
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
legend('Raw Y Data', 'Fitted Y Data'); %Add Legend for the two lines
title('Chi Squared with fitted curve');
text = 'Chi Squared: ' + string(ChiSquared); % make string for textbox
annotation('textbox',[.2 .9 0 0],'String',text,'FitBoxToText','on'); % annotate graph with textbox
