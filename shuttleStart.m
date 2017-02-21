

%shuttleStart is a script to initialise the simulation of the heat flow
%through a space shuttle tile.

%Construct the Input boxes and default values for the various tile
%simulation parameters

prompt          = {'Simulation Maximum Time (s):'...
                   'Number of Time Steps:'...
                   'Simulation Tile Thickness (m):'...
                   'Number of Spacial Steps:'...
                   'Simulation Method (Forwards/Backwards/Dufort/Crank)'...
                   'Space Shuttle Data Point'...
                   'Graph? (True/False)'};
name            = 'Insert Tile Simulation Variables';
defaultAnswer   = {'4000','501','0.05','21','crank','597','true'};
numlines        = 1;

%Bring up picture interface so that people can see the tile selection
scrsz = get(groot,'ScreenSize');        %Finds out the size of the current screen
figure('Position',[1 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/2]) %Aligns Image to the top right of screen
img =imread('fullShuttleMerged.jpg');
image(img)

%Register the inputs from the input box
inputs = inputdlg(prompt,name,numlines,defaultAnswer);
close %Close Picture
tMax    =str2double(inputs{1});
nt      =str2double(inputs{2});
xMax    =str2double(inputs{3});
nx      =str2double(inputs{4});
method  =(inputs{5});
point   =str2double(inputs{6});
doPlot  =(inputs{7});

%Begin Basic Error Checking
if tMax < 0
    error('Please enter a positive simulation time input')
    
elseif nt > 4000
    error('Lots of Timesteps, Simulation will be slow')
    
elseif xMax >1
    error('More than 1m thick tiles? I''surprised the rocket took off in the first place')
    
end

%Call the Shuttle Simulation Function 
shuttleSimulation1D(tMax, nt, xMax, nx, method, doPlot, point)


