function [x t u] = shuttleSimulation1D(tmax, nt, xmax, nx, method, doplot,point)
% Function for modelling temperature distribution in 1D in a space shuttle tile
% D N Johnston  10/04/2016
% Amended by N Rupp 2016
%
% Input arguments:
% tmax   - maximum time
% nt     - number of timesteps in simulation
% xmax   - total thickness of space shuttle tile
% nx     - number of spatial steps in simulation
% method - chosen simulation solution method ('forward', 'backward' etc)
% doplot - true to plot graph; false to suppress graph.
% point  - point on the space shuttle (tile No)
%
% Return arguments:
% x      - distance vector
% t      - time vector
% u      - temperature matrix
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = shuttle(4000, 501, 0.05, 21, 'forward', true, 597);

%Defines the space shuttle tile's  properties
thermcon    = 0.142;   % Thermal Conductivity W/(m K)
density     = 352;     % Density converted from 22 lb/ft^3 to (kg/m^3)
specheat    = 1256;    % 0.3 Btu/lb/F at 500F

%Calculates the value for alpha used in the Fourier equation
alpha = thermcon /(density * specheat);

%Load the tile property's raw data from an analysed Data file
fileName= sprintf('%d.mat',point');
load (fileName)

%CALCULATE VALUES FROM RAW DATA
timeMax = tmax;          %Extracts the total time of the raw data
dt = timeMax / (nt-1);  % Calculates the timestep dt
t  = (0:nt-1)*dt;       %generates a time vector
dx = xmax/(nx-1);       %Calculates the distance step by dividing total distance by no steps
x  = (0:nx-1)*dx;       %generates a distance vector
u  = zeros(nt, nx);     %Generates an empty matrix which looks at each time and temp step
p  = alpha * dt / dx^2; %Calculates the p value
ivec = 2:nx-1;          % Creates a vector of indices for internal points

%Set initial conditions to 60F throughout (converted to celsius)
u([1 2], :) = (60-32)*5/9;

% Begin the Main timestepping loop to model the heat flow through a tile
% using various simulation methods
for n = 2:nt - 1
    if t(n)<2000
      % Linearly interplolate outer surface data where required
        tempf = interp1(timeData, tempData, t(n+1), 'linear', 'extrap');
    else
        % To prevent interpolation after the end of the data set the tile
        % exterior to 60F
        tempf = 60;
    end
        % Convert from degrees F to degrees C
        R = (tempf - 32) * 5 / 9;
        % Set as outer surface RHS boundary condition
        u(n+1, nx) = R;
        
    % Select method for simulation.
    switch method
        case 'forwards'
            
           
            % LHS Neumann boundary condition to represent zero heat flow in
            % the interior of the rocket
            u(n+1, 1) = (1-2*p)*u(n,1)+ 2*p*u(n,2);
            % Internal points
            ivec=2:nx-1;
            u(n+1, ivec) = (1 - 2 * p) * u(n, ivec) + ...
                p * (u(n, ivec-1) + u(n, ivec+1));
           
            % RHS Dirichlet boundary condition of the exterior surface
            u(n+1, nx) = R;
            
            
        case 'dufort'
            % LHS Neumann boundary condition to represent zero heat flow in
            % the interior of the rocket 
            u(n+1, 1) = ((1-2*p)*u(n-1,1)+ 4*p*u(n,2))/(1+2*p);
            % Internal points
            ivec=2:nx-1;
            u(n+1, ivec) = ((1 - 2 * p) * u(n-1, ivec) + ...
                2* p * (u(n, ivec-1) + u(n, ivec+1)))/(1+2*p);
            % RHS Dirichlet boundary condition
            u(n+1, nx) = R;
            
        case 'backwards'
            %Backwards Differencing Method
            
            %Initalise the Initial Vlues of the backwards differencing
            %matrix
            % [b c       | [u|   [d|
            % |a b c     | |u|   |d|
            % |   - - -  | |-| = |-|
            % |     a b c] |u]   |d]
            
            %Initial Values of b,c,d in first row
            b(1)      = 1+2*p;
            c(1)      = -2*p;
            d(1)      = u(n,1); %LHS Neumann Boundary Condition
            
            %General Formulae for subsequent values of
            a(2:nx-1) = -p;
            b(2:nx-1) = 1 + 2*p;
            c(2:nx-1) = -p;
            d(2:nx-1) = u(n, 2:nx-1);
            
            %Final and fixed values for a,b,d in the last row
            a(nx)     = 0;
            b(nx)     = 1;
            d(nx)     = R; %RHS Boundary Condition
            
            %Tri-Diagonal Matrix Method is then called to solve the matirx for the
            %temperature vector u (located at the bottom of the code)
            u(n+1,:) = tdm(a,b,c,d);
            
        case 'crank'
            %Crank Nicholson Method
            %Approximation based on Fourier's equation at the mid-point
            %between two timesteps, (x_i, t_n+delta_t/2)
            
            %Second order accuracy in time and space, therefore more
            %accurate but more complicated then backwards differencing,
            %therefore more computation time.
            
            %Initalise the Initial Values of the backwards differencing
            %matrix in the form:
            % [b c       | [u|   [d|
            % | a b c    | |u|   |d|
            % |   - - -  | |-| = |-|
            % |     a b c] |u]   |d]
            
            %Initial Values of b,c,d in first row
            b(1)      = 1+p;
            c(1)      = -p;
            d(1)      = (1-p)*u(n,1) + p*u(n,2); %LHS Neumann Boundary Condition
            
            %General Formulae for subsequant values of
            a(2:nx-1) = -p/2;
            b(2:nx-1) = 1 + p;
            c(2:nx-1) = -p/2;
            d(2:nx-1) =(p/2)*u(n,1:nx-2)+(1-p)*u(n,2:nx-1)+(p/2)*u(n,3:nx);
            
            %Final and fixed values for a,b,d in the last row
            a(nx)     = 0;
            b(nx)     = 1;
            d(nx)     = R; %RHS boundary condition
            
            %Tri-Diagonal Matrix Method is then called to solve the matrix for the
            %temperature vector u (Located at bottom of code)
            u(n+1,:) = tdm(a,b,c,d);
            
        otherwise
            %Error catching on incorrect method entries
            error (['Undefined method: ' method])
            return
    end
end

%Begin Plotting
if doplot
    
    %Find out pixel size of the screen
    scrsz = get(groot,'ScreenSize');
    %Initialise a figure in the top left of the user's screen
    figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/2])
    
    %Generate a surface plot of x,t,u from the simulation
    surf(x, t, u); 
    shading interp %Shading parameters
    
    %Switch statement to extract Method name for graph title
    switch method
        case 'forwards'
            method= 'Forwards Differencing';
        case 'dufort'
            method= 'Dufort-Frankel';
        case 'backwards'
            method= 'Backwards Differencing';
        case 'crank'
            method= 'Crank-Nicolson';
    end
    
    %Generate graph title from simulation variables using sprintf
    graphTitle = sprintf('1D heat flow through insulated shuttle tile no.%d\n Method: %s, Max Time: %ds, Thickness: %gm\nTime Step: %gs, Spacial Step: %gm',point,method,tmax,xmax,dt,dx);
    title(graphTitle);
    
    %Label the Axes
    xlabel('Distance from Inner Surface, x (m)');
    ylabel('Time, t (s)');
    zlabel('Tile Temperature, T (^\circC)');
    
    
    %Open up a message box that, when clicked, closes the graph and returns
    %to main menu
    uiwait(msgbox('Click to Return to Parameter Menu'))
    close
    shuttleStart
end
% End of shuttle function


% =========================================================================
% Tri-diagonal matrix solution
function x = tdm(a,b,c,d)
n = length(d);

% Eliminate a from the matrix terms
for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

%Solve the first term in the matrix
x(n) = d(n) / b(n);

% Loop backwards to find other x values by back-substitution
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end
% =========================================================================

