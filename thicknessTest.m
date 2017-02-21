function minThickness = thicknessTest(startThickness,endThickness,step)
%The function thickness test investigates the effect of thickness on the
%temperature of the internal surface of an insulated heat shield.


%Set Up Parameters of the Prompt Box
prompt          = {'Start Thickness (m):'...
                   'End Thickness (m):'...
                   'Thickness Size Increment (m):'};              
name            = 'Thickness Testing';
defaultAnswer   = {'0.01','0.13','0.02'};
numlines        = 1;

%Call the prompt box
inputs = inputdlg(prompt,name,numlines,defaultAnswer);

%Initialise the parameters of the thickness test with the answers from the
%input dialogue.
startThickness  =str2double(inputs{1}); %Smallest test value of thickness
endThickness    =str2double(inputs{2}); %Largest test value of thickness
step            =str2double(inputs{3}); %test increment value

%Initialise constant variables to be used in the 1D Shuttle Simulation
%function during the testing
tMax    = 4000;       %maximum time
nt      = 501;        %number of timesteps in simulation
nx      = 21;         %number of spatial steps in simulation
method  = 'crank';    %chosen simulation solution method ('forward', 'backward' etc)
doPlot  = false;      %true to plot graph; false to suppress graph.
point   = 597;        %Space shuttle tile point chosen


%Calculate graph parameters for plotting purposes
dt = tMax / (nt-1);
noSteps= floor((endThickness-startThickness)/step);

%Initialise basic graph parameters for looped plotting
figure(1);
xlim([0 4000])

%Intialise Iterator for vector indexing purposes
i=1;        
hold off % Prevent old graphs remaining on the figure

for testThickness= startThickness:step:endThickness
    hold on %Retain Graph eatch loop
    
    %Call Shuttle function for the current test thickness
    [~, t, u] = shuttleSimulation1D(tMax, nt, testThickness, 21, 'crank', false, 597);
    
    %Plot the tile interior's temperature distribution and save the current
    %thickness to a variable h for the legend.
    h(i) = plot (t,u(:,1),'DisplayName',sprintf('%gm',testThickness));
    
    %Prevents unneccessary calculations if the temperature no longer has an
    %effect on the interior temperature.
    if u(1,1)==u(end,1)
        minThickness=testThickness;
        break
    end
    i=i+1;
end

%Plot the original temperature on the tile's exterior as a comparison point
outerU=u(:,end);
h(end+1)= plot(t,outerU,'k--','DisplayName','Outer Temp');

%Prepare and initialise the graph title
graphTitle = sprintf('Insulated Temperature Variation at Different Thicknesses in Tile no.%d\n Method: %s, Max Time: %ds, Time Step: %gs',point,method,tMax,dt);
title(graphTitle);

%Initialise Legend and scale and label axes
legend(h)
xlabel('Time, t (s)')
ylabel('Temperature, T (^\circC)')
ylim([0 1000]);

%Prevent future graphs from being affected by this graph
hold off

end

