% 1D 1 clast conduction only
% Clear all figures 
clf;
clear all;
 
clastnum = 1; % Number of clasts in spatter pile
 
% Define model parameters
% Naming scheme: x is a location. dx is a thickness. 
% Node is nodal location (when combined with x) and number of nodes (when
% combined with dx)
% Should only change dxclast out of all model parameters 
dysystem = 4.; % (m) Thickness of entire system (area shown in model window)
dysystemnode = (801); % Number of nodes for the thickness of the entire system. Controls resolution.
dystep = dysystem/(dysystemnode-1); % (m) Length one grid step (dx)
dyclast = 0.05; % (m) Thickness of one clast
dyclastnode = dyclast/dystep; % Number of nodes for thickness of one clast
ypilebottom = 1.; % (m) Boundary between the ground and the bottom of the spatter pile
ypilebottomnode = ((ypilebottom/dystep)+1); % Node of the boundary between the ground and the bottom of the spatter pile
yclastonetop = (ypilebottom+dyclast); % (m) Boundary between the top of the clast and the air (if spatter pile has only one clast)
yclastonetopnode = (ypilebottomnode+dyclastnode); % Node of the boundary between the top of the clast and the air (if spatter pile has only one clast)
yclasttwotop = (yclastonetop + dyclast);
yclasttwotopnode = (yclastonetopnode + dyclastnode); 
yclastthreetop = (yclasttwotop + dyclast);
yclastthreetopnode = (yclasttwotopnode + dyclastnode);
dypile = (dyclast*clastnum); % (m) Thickness of spatter pile
dypilenode = (dyclastnode*clastnum); % Number of nodes for thickness of spatter pile
jstep = 1; % step after every flow
 
% Thermal constants
kground = 1.7082; % (W/(mK)) Thermal conductivity of the ground (previously deposited basalt)
kclast = 1.7082; % (W/(mK)) Thermal conductivity of the clast (basalt)
kair = 0.0264; % (W/(mK)) Thermal conductivity of the air
cpground = 1490.; % (J/(kgK)) Specific heat
cpclast = 1490.;
cpair = 1006.;
rhoground = 2910.; % (kg/m3) Density 
rhoclast = 2910.;
rhoair = 1.204; 
rhocpground = rhoground*cpground; % (J/(K(m^3))) Volumetric heat capacity
rhocpclast = rhoclast*cpclast;
rhocpair = rhoair*cpair;
h = 8.; % (W/m^2 K) Heat transfer coefficient
kappaground = kground/rhocpground; % (m^2/s) Thermal diffusivity
kappaclast = kclast/rhocpclast; % (m^2/s) Thermal diffusivity
kappaair = kair/rhocpair; % (m^2/s) Thermal diffusivity
emis = 0.96; % emissivity (unitless)
sb = 0.000000056703; % W/m^2K^4 stefan-boltzmann constant
 
% Timestep
dtexp = dystep^2/2/kappaclast; % (s) Limitation for explicit timestep 
dt = dtexp./10; % (s) Length of one timestep
 
% Things to vary with different model runs
 
tnum = 50000; % Number of timesteps (How many timesteps do you want to run the code for?)
%tclasttwo = 10.; % (s) Time when a new clast is deposited
%tclastthree = 20.; 
%tbetween = 100.; % (s) Time between subsequent clasts deposited
 
% Setting up the vector for nodal point positions (the righthand size of our matrix system is the vector) 
y=0:dystep:dysystem; % The vector is essentially one grid step by the thickness of the system (0:0.005:2)
 
% Define the initial temperature profile for the system
Tground = 293.; % (K) Initial ground temperature
Tclast = 1373.; % (K) Initial clast temperature
Tair = 293.; % (K) Initial temperature of the air
Tbackground = 293.; % (K) Initial temperature background for the entire system
 
% Tell the model to set the initial temperature profile for the system
for i=1:1:dysystemnode % For all locations of i in the system
    Ti(i) = Tbackground; % We say that initially the entire system is at 293K
if((i)>= ypilebottomnode+1 && (i)<= yclastonetopnode)
    Ti(i) = Tclast; % If we are located within the clast (including bottom node and top node of the clast) our initial temp is now 1403K
end 
end 
 
% Time cycle
timesum=0; % How much time has elapsed (we're still looking at what the model is initially, so no time has yet passed)
 
% Now let's populate the matrix
L = sparse(dysystemnode,dysystemnode); % The lefthand side of the matrix is 401 nodes by 401 nodes and is mostly made up of zeros (sparse)
R = zeros(dysystemnode,1); % The righthand side of the matrix is the vector and it is 401 nodes by 1
 
b=1;
for t=1:1:tnum
    
    % Bottom of system held at constant temp - Tground
    L(1,1) = 1; % First point in the matrix is 1
    R(1,1) = Tground; % First point in the vector is Tground
 
    
    if jstep == 1 % If only one clast in the system
        
        % Define internal nodes of the ground (conduction only)
    for i=2:1:ypilebottomnode
        L(i,i-1) = -kclast/(rhoclast*cpclast*(dystep^2));
        L(i,i) = (kclast/(rhoclast*cpclast*(dystep^2))) + (kclast/(rhoclast*cpclast*(dystep^2))) + 1/dt;
        L(i,i+1) = -kclast/(rhoclast*cpclast*(dystep^2));
        R(i,1) = Ti(i)/dt;
    end
    
        % Define ground/clast boundary (conduction only) (part of the
        % clast)
  
        L(ypilebottomnode+1,ypilebottomnode) = -kclast/(rhoclast*cpclast*(dystep^2));
        L(ypilebottomnode+1,ypilebottomnode+1) = (kclast/(rhoclast*cpclast*(dystep^2))) + (kclast/(rhoclast*cpclast*(dystep^2))) + 1/dt;
        L(ypilebottomnode+1,ypilebottomnode+2) = -kclast/(rhoclast*cpclast*(dystep^2));
        R(ypilebottomnode+1,1) = Ti(ypilebottomnode+1)/dt;
        
   
        
        % Define internal nodes of the clast (conduction only)
    for i=ypilebottomnode+2:1:yclastonetopnode-1
        L(i,i-1) = -kclast/(rhoclast*cpclast*(dystep^2));
        L(i,i) = (kclast/(rhoclast*cpclast*(dystep^2))) + (kclast/(rhoclast*cpclast*(dystep^2))) + 1/dt;
        L(i,i+1) = -kclast/(rhoclast*cpclast*(dystep^2));
        R(i,1) = Ti(i)/dt;
    end
    
        % Define clast/air boundary (conduction, convection, radiation)
    
        L(yclastonetopnode,yclastonetopnode-1) = -kclast/(rhoclast*cpclast*(dystep^2));
        L(yclastonetopnode,yclastonetopnode) = (kclast/(rhoclast*cpclast*(dystep^2))) + (kair/(rhoclast*cpclast*(dystep^2))) + 1/dt;
        L(yclastonetopnode,yclastonetopnode+1) = -kair/(rhoclast*cpclast*(dystep^2));
        R(yclastonetopnode,1) = Ti(yclastonetopnode)/dt;
            
        
        % Define internal nodes of the air (conduction only)
    for i=yclastonetopnode+1:1:dysystemnode-1
        L(i,i-1) = -kair/(rhoair*cpair*(dystep^2));
        L(i,i) = (kair/(rhoair*cpair*(dystep^2))) + (kair/(rhoair*cpair*(dystep^2))) + 1/dt;
        L(i,i+1) = -kair/(rhoair*cpair*(dystep^2));
        R(i,1) = Ti(i)/dt;
    end 
    
    % Top of system held at constant temp - Tair
   
    L(dysystemnode,dysystemnode) = 1; % Last point in the matrix is dysystemnode
    R(dysystemnode,1) = Tair; % Last point in the vector is Tair
    end
     
% Solution
Tnew = L\R;
 
timevec(b)=timesum;
tempvec(b)=Ti(yclastonetopnode-1);
b=b+1;
 
disp(timesum);
disp(jstep);
 
% Open figure
% Plotting implicit solution 
subplot(2,1,1);
plot(Ti,y,'k');
axis([0.9*Tground 1.1*Tclast 0 dysystem]);
title(['Implicit: step=',num2str(t),' time,sec=',num2str(timesum)]);
ylabel('y, m');
xlabel('Temperature, K');
% Plot cooling rate
subplot (2,1,2);
plot (timesum, Ti(yclastonetopnode-1),'o'); %This is the "location" of the thermocouple, we want it to be one node below the top of the first clast
hold on;
axis([0 500 300 1410]);
xlabel('Time, sec');
ylabel('Temperature, K');
text (10,10,'1025 sec - 2nd flow');
 
% Stop for 0.1 second
drawnow % draw above figure without delay
pause(0.1);
 
% Stop for 0.1 second
drawnow % draw above figure without delay
pause(0.1);
 
% Time counter
timesum = timesum+dt;
 
 
% Reassign temp profile for next round
Ti = Tnew;
 
end
 
    b=1;
for i=1:length(timevec)
    if tempvec(i)>=900
        Time(b)=timevec(i);
        Temp(b)=tempvec(i);
    end
    b=b+1;
end
 
format longg
 
TimeFinal=Time(end)
TempFinal=Temp(end)
