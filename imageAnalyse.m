function usableData = imageAnalyse(fileName, fileType) 
%Image Analyse Function
%Imports the image that is to be analysed.
%This function is an image reading function that was designed to convert a
%.jpg picture of a graph with red data points into scaled data points for use in MATLAB
%This version is specifically adapted to read .jpg
%
%It:
%   1) Imports image and gets information about the size of the JPEG
%   2) Extracts the relative location of each red data point, the graphs
%   origin and axis length
%   3) Scales down the relative lengths to give correctly scaled values of temperature
%   and time 
%   4) Stores the extracted data in a .m file.
%   5) Displays Graph of imported data 
%
%Inputs: 
%   fileName: The name of the graph's image.
%   fileType: The type of file.
%Outputs:
%   .m file containing the extracted data for use in other functions

%DUMMY VALUES FOR TESTING
%fileName='temp597';
%fileType='.jpg';
%Set Up Parameters of the Prompt Box
prompt          = {'File Name:'...
                   'File Type:'};
                   
name            = 'Image Analyser';
defaultAnswer   = {'597','.jpg'};
numlines        = 1;

%Call the prompt box
inputs = inputdlg(prompt,name,numlines,defaultAnswer);

%Initialise the parameters of the thickness test with the answers from the
%input dialogue.
fileName  =inputs{1}; %Smallest test value of thickness
fileType  =inputs{2}; %Largest test value of thickness

%Extract all the image information from the function inputs
img=imread([fileName fileType]);

%Analyse the image to extract the image dimensions
imageSize   =size(img);
numberRows  =imageSize(1);
numberCols  =imageSize(2);

%Find the Coordinates of the Red Data Points
    %RGB Values for the target colours for the red data points
    r=150;
    g=120;
    b=120;

    %Begin a for loop to find all the red data points
    %It loops through every pixel column looking for red points. When it finds
    %multiple red point it takes an average of the point's position and stores
    %the information in rawGraphData.
    for i=1:numberCols
        rawGraphData(i)=mean(find(img(:,i,1)>r & img(:,i,2)<g & img(:,i,3)<b));
    end

%Find the Location of the Axes
    %RGB Values for the black axis points search critera
    r=200;
    g=100;
    b=100;

    %Initialise a zeros matrix for data information to be inputted into
    y= zeros(numberRows,numberCols);

    %Begin the for loop to find the location of all black pixels in the picture
    %It generates a matrix the same dimension as the imported picture. Each
    %Each pixel in  the image is then assigned a value:
    %as 1 (black) or 0 (Not black)
    %One limitation is that this for loop assumes a landscape graph 
    for i=1:numberCols
        if i<=numberRows
            y(:,i)=(img(:,i,1)<r & img(:,i,2)<g & img(:,i,3)<b);
            ySumVertical(:,i)=sum(y(:,i));
            ySumHorizontal(i,:)=sum(y(i,:));
        else
            y(:,i)=(img(:,i,1)<r & img(:,i,2)<g & img(:,i,3)<b);
            ySumVertical(:,i)=sum(y(:,i));

        end
    end

    %Find coordinates of the Origin
    %Looks for the row and column where there are the most black pixels in a
    %line in the picture.
    [xBlackPoints originX] = max(ySumVertical(:));
    [yBlackPoints originY] = max(ySumHorizontal(:));

%Find the extent of the x axis
    %Begin the iteration on the right y coordinate from the origin

    %Initialise values for the while loop
    i=originX;
    test=1;
    %The while loop searches for adjacent black pixels. It begins at the axis
    %origin and searches out until there are no longer black pixels, At which point it knows the length of the axis in pixels. 
    while test
        if i==numberCols
            %An error catching loop to prevent it looping past the extent of
            %the picture
            break
        else
            %Updates the test variable. if the pixel is black Test== 1,
            %non-black == 0 
            test=y(originY,i);
            %Add one to the axis length
            i=i+1;
        end
    end
    %Final Value for the X-Axis Length
    xAxisEnd= i;

%Generate the length of the Y Axis in a similar manner to the X-Axis
    i=originY;
    test=1;
    %Loop from the origin
    while test
        if i==0
            break
        else
    test=y(i,originX);
    i=i-1;
        end
    end
    yAxisEnd= i;

axisData=[originX originY xAxisEnd yAxisEnd];

%Plot raw Graph data to check for consistency and axis scaling
    amendedData=amendData(rawGraphData,axisData);
    graphTitle= sprintf('Raw Data for Tile Number %s',fileName);
    title(graphTitle)
    hold on
    plot(amendedData(1,:),amendedData(2,:));
    xlabel('Time (s)');
    ylabel('Temperature (F)');
    timeData=amendedData(1,:);
    tempData=amendedData(2,:);

%save data to .mat file with same name as image file
save(fileName, 'timeData', 'tempData')

end





%=========================================================================
function amendedData = amendData(rawData,axisData)

%Extract Variables from function input to make code more readable
originX=axisData(1);
originY=axisData(2);
maxTimeCoord=axisData(3);
maxTempCoord=axisData(4);


%Readings from the maximum increments of the graphs
maxTimeValue=2000; % Maximum time axis value in seconds (s)
maxTempValue=2000; % Maximum temp axis value in farenheight (F)


%Find out the absolute pixel lengths of the time and temp axes
timeAxisLength=abs(maxTimeCoord-originX);
tempAxisLength=abs(maxTempCoord-originY);

%Calculate the Scale Factors between pixel length and actual length
timeScaleFactor=maxTimeValue/timeAxisLength;
tempScaleFactor=maxTempValue/tempAxisLength;

%Create and scale time vector
timeVector=[0:timeAxisLength];
timeAmend=timeScaleFactor.*timeVector;

%Create and scale temp vector 
tempVector=rawData(originX:(originX+timeAxisLength));
tempVector=originY-tempVector;
tempAmend=tempScaleFactor.*tempVector;

%Clears out the NaN values from the temperature vector to prevent
%interpolation issues later on.
tempLength=length(tempAmend);
for n= tempLength:-1:1
    if isnan(tempAmend(n))
       tempAmend(n)=[];
       timeAmend(n)=[];
    end
end

% Puts the final time and temp axis data into one matrix for outputting
% from the function
amendedData(1,:)=timeAmend;
amendedData(2,:)=tempAmend;
end

