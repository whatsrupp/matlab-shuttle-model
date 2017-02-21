function [designThickness] =thicknessShooter (targetTemp)
%Thickness Shooter is a function that calculates the minimum allowable
%thickness of the space shuttle given the boundary condition of maximum
%allowable temperature

%Initialise prompt box parameters
prompt          = {'What is the Maximum Allowable Temperature on the Inner Surface (degrees C)?'}          
name            = 'thicknessShooter';
defaultAnswer   = {'200'};
numlines        = 1;

%Call the prompt box
inputs = inputdlg(prompt,name,numlines,defaultAnswer);

%Initialise the parameters of the thickness test with the answers from the
%input dialogue.
targetTemp =str2double(inputs{1}); %Maximum temperature boundary condition

%Fixed Shuttle Simulation Parameters
tMax= 4000;         %maximum time
nt = 501;           %number of timesteps in simulation
nx = 50;            %number of spatial steps in simulation
method= 'crank';    %chosen simulation solution method ('forward', 'backward' etc)
doPlot= false;      %true to plot graph; false to suppress graph.
point   = 597;      %Space shuttle tile point choice

%Initial Thickness Guesses
guess1=0.1;
guess2=0.01;

%Begin Shooting Method Loop
for i=1:50
    
    %Call Shuttle Simulation to get the temperature matrices of each
    %simulation
    [~,~,u1] = shuttleSimulation1D(tMax, nt, guess1, nx, method, doPlot, point);
    [~,t,u2] = shuttleSimulation1D(tMax, nt, guess2, nx, method, doPlot, point);
    
    %Find the temperature of the insulated interior
    maxTemp(1)= max(u1(:,1));
    maxTemp(2)= max(u2(:,1));
    
    %Calculate the error for the current loop
    error = targetTemp-maxTemp(2);
    
    %Test for error each loop, if the criteria are met then break the loop
    if abs(error)<1
        break
    end
    
    %Call Shooting Method 
    guess3 = guess2 + error*((guess2-guess1)/(maxTemp(2)-maxTemp(1)));
    
    %Update guesses for next guessing loop
    guess1 = guess2;      
    guess2 = abs(guess3);
    
end

%State result from the shooting method
minThickness = guess2;

%Begin Plotting the Graph
%Variable h is used to store different legend entries
hold off
h(1)=plot(t,u2(:,1),'DisplayName','Inner Temp Variation');
hold on
h(2)=plot(t,u2(:,end),'DisplayName', 'Outer Temp Variation');

%Create a comparison line between the shooter and the final thickness value
testSize=size(t);
testLine(1:testSize(2))=targetTemp;
h(3)=plot(t,testLine,'k--','DisplayName','Target Temp Line');

%Initialise and Call Graph Title
graphTitle= sprintf('Inner Temperature Shooting Method\n Target Temperature: %dC,  Minimum Thickness: %gm',targetTemp,minThickness);
title(graphTitle)

%Initialise Legend and Axis Labelling
legend (h);
ylabel ('Temperature ( ^\circC)');
xlabel ('Time (s)');
