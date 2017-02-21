% Script to plot image of measured temperature, and trace it using the mouse.
%
% Image from http://www.columbiassacrifice.com/techdocs/techreprts/AIAA_2001-0352.pdf
% Now available at http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.26.1075&rep=rep1&type=pdf

%Imports jpg image data in variable img
name = '597';
img=imread([name '.jpg']);

%Creates Completely Blank Figure Window
figure (4);
%Plots image onto the figure
image(img);
%Keeps the image on the screen
hold on

%Uses Title As the First User Command
title('Left Click Curve Profile, right Click to end data collecting')

%Initialise Time, Temperature and Axis Information Vectors
timeData = [];
tempData = [];
timeAxisData = [];
tempAxisData = [];

%Commence infinite loop for the manual entry of data into clicking the
%graph
while 1 % infinite loop
    
    % get one point using mouse
    [x, y, button] = ginput(1);
    
     % break if anything except the left mouse button is pressed
    if button ~= 1
        break
    end
    %Plot a green circle when the left mouse button is clicked
    plot(x, y, 'og') % 'og' means plot a green circle.
    
    % Add data points to vectors. 
    timeData = [timeData, x];
    tempData = [tempData, y];
end

%Collect Data Calibration Data

title('Click on Origin');
[x, y, button] = ginput(1);
plot(x, y, 'og');
timeAxisData = [timeAxisData,x];
tempAxisData = [tempAxisData,y];

%Collect Temperature Scaling Factor
title('Click on 1000F Temp axis Point');
[x, y, button] = ginput(1);
plot(x, y, 'og');
tempAxisData = [tempAxisData,y];

%Collect Time Scaling Factor
title('Click on 500s Time axis Point');
[x, y, button] = ginput(1);
plot(x, y, 'og');
timeAxisData = [timeAxisData,x];

%Initialise Alignment variables for code clarity
timeAlign = timeAxisData(1);
tempAlign = tempAxisData(1);

%Initialise scaling variables for code clarity
timescale= 500/(timeAxisData(2)-timeAxisData(1));
tempscale= 1000/abs((tempAxisData(2)-tempAxisData(1)));

%Realign and then Rescale the Temperature data to give correct graph values
%in farenheight
tempData= abs(tempData-tempAlign);
tempData= tempscale.*tempData;

%Realign and then Rescale the Time data to give correct graph values in
%seconds
timeData= timeData-timeAlign;
timeData= timescale.*timeData;

% sort data and remove duplicate points.
[timeData, index] = unique(timeData);
tempData = tempData(index);

%save data to .mat file with same name as image file
save(name, 'timeData', 'tempData')
hold off

