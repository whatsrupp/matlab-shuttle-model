
%Stability test is a script that runs through the shuttle simulation
%numerous times using all 4 simulation methods. Each iteration the timestep
%changes. The output graph gives insight into the effectiveness of each
%simulation as the inputted number of timesteps change.

prompt          = {'Minimum Number of Timesteps:'...
                   'Maximum Number of Timesteps:'...
                   'Timestep Number Increment (s):'};
                   
name            = 'Stability Testing';
defaultAnswer   = {'41','1001','50'};
numlines        = 1;

%Call the prompt box
inputs = inputdlg(prompt,name,numlines,defaultAnswer);

%Initialise the parameters of the thickness test with the answers from the
%input dialogue.
minTimeStep  =str2double(inputs{1}); %Minimum amount of time steps
maxTimeStep  =str2double(inputs{2}); %Maximum amount of time steps
increment    =str2double(inputs{3}); %test increment value

%Set the Shuttle simulation parameters
tMax= 4000;         %maximum time
nt = 501;           %number of timesteps in simulation
xMax=0.05;          %Overall thickness of the tile
nx = 21;            %number of spatial steps in simulation
doPlot= false;      %true to plot graph; false to suppress graph.
point   = 597;      %Space shuttle tile point chosen


%initialise the programme with i=0 at the outer boundary
i=0;

%For a an increasing number of timesteps
for nt = minTimeStep:increment:maxTimeStep 
    %Increment to the next step in the displacement matrix
    i=i+1; 
    
    %Time step at point 1= maximum time/ number of steps -1
    dt(i) = tMax/(nt-1); 
    %Displays the value of the things for increasing time step
    disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s']) 
    [~, ~, u] = shuttleSimulation1D(tMax, nt, xMax, nx, 'forwards', false, point); 
    uf(i) = u(end, 1); 
    [~, ~, u] = shuttleSimulation1D(tMax, nt, xMax, nx, 'backwards', false, point);
    ub(i) = u(end, 1); 
    [~, ~, u] = shuttleSimulation1D(tMax, nt, xMax, nx, 'crank', false, point);
    uc(i) = u(end, 1); 
    [~, ~, u] = shuttleSimulation1D(tMax, nt, xMax, nx, 'dufort', false, point);
    ud(i) = u(end, 1); 
end 

%Generate the data for the real temperature line on the graph to act as a
%comparison point.
testSize = size (ud);
testConstant = uc(1,3);
testLine(1:testSize(2))= testConstant;

%Plot Graph 
hold off %Clear previous graph
 
plot(dt, [uf; ub; uc; ud]); %Plot the test results
hold on                     %Prevent graph clearing
plot(dt, testLine, 'k--');  %Plot the real temperature value


%initialise and call the graph title using sprintf
graphTitle=sprintf('Stability Graph\n Min No Timesteps = %d, Max No Timesteps = %d, Increment = %d',minTimeStep,maxTimeStep,increment);
title(graphTitle)

%Label and Scale Axes
xlim('auto')
ylim([140 180])
xlabel('Size of Time Steps (s)');
ylabel('Final Temp (^\circC)');

%Locate Legend in the Top Middle
legend ('Forward','Backward','Crank-Nicolson','Dufort Frankel','Actual Temperature','location','north')
hold off