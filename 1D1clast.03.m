% 1D multiple clasts CCR
% hightemplowv

% Clear all figures 
clf;
clear all;

clastnum = 1; % Number of clasts in spatter pile
dyclast = 0.03; % (m) Thickness of one clast

tnum = 30000; % Number of timesteps (How many timesteps do you want to run the code for?)
tclasttwo = 5.; % (s) Time when a new clast is deposited
tclastthree = 10.; 


% Define model parameters
% Naming scheme: x is a location. dx is a thickness. 
% Node is nodal location (when combined with x) and number of nodes (when
% combined with dx)
% Should only change dxclast out of all model parameters 
dysystem = 4.; % (m) Thickness of entire system (area shown in model window)
dysystemnode = (801); % Number of nodes for the thickness of the entire system. Controls resolution.
dystep = dysystem/(dysystemnode-1); % (m) Length one grid step (dx)
dyclastnode = dyclast/dystep; % Number of nodes for thickness of one clast
ypilebottom = 1.; % (m) Boundary between the ground and the bottom of the spatter pile
ypilebottomnode = ((ypilebottom/dystep)+1); % Node of the boundary between the ground and the bottom of the spatter pile
yclastonetop = (ypilebottom+dyclast); % (m) Boundary between the top of the clast and the air (if spatter pile has only one clast)
yclastonetopnode = (ypilebottomnode+dyclastnode); % Node of the boundary between the top of the clast and the air (if spatter pile has only one clast)
yclasttwotop = (yclastonetop + dyclast);
yclasttwotopnode = (yclastonetopnode + dyclastnode); 
yclastthreetop = (yclasttwotop + dyclast);
yclastthreetopnode = (yclasttwotopnode + dyclastnode);
%xpiletop = (xpilebottom+(dxclast)); % (m) Boundary betweem the top of the spatter pile and the air at the start
%xpiletopnode = (xpilebottomnode+(dxclastnode)); % Node of the boundary between teh top of the spatter pile and the air at the start
%xpiletopnode = 259;
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



% Setting up the vector for nodal point positions (the righthand size of our matrix system is the vector) 
y=0:dystep:dysystem; % The vector is essentially one grid step by the thickness of the system (0:0.005:2)

% Define the initial temperature profile for the system
Tground = 293.; % (K) Initial ground temperature
Tclast = 1373.; % (K) Initial clast temperature
Tair = 293.; % (K) Initial temperature of the air
Tbackground = 293.; % (K) Initial temperature background for the entire system

% Tell the model to set the initial temperature profile for the system
for i=1:1:dysystemnode % For all locations of i in the system
    T0imp(i) = Tbackground; % We say that initially the entire system is at 293K
if(y(i)> ypilebottom && y(i)<= yclastonetop)
    T0imp(i) = Tclast; % If we are located within the clast (including bottom node and top node of the clast) our initial temp is now 1403K
end 
end 



% Time cycle
timesum=0; % How much time has elapsed (we're still looking at what the model is initially, so no time has yet passed)

% Now let's populate the matrix
L = sparse(dysystemnode,dysystemnode); % The lefthand side of the matrix is 401 nodes by 401 nodes and is mostly made up of zeros (sparse)
R = zeros(dysystemnode,1); % The righthand side of the matrix is the vector and it is 401 nodes by 1

% Define the first point in the matrix and the first point in the vector
L(1,1) = 1; % First point in the matrix is 1
R(1,1) = Tground; % First point in the vector is 293K

% Tells the model to add a new clast to the spatter pile
b=1;
for t=1:1:tnum
    %xpiletopnode = (xpilebottomnode+(dxclastnode*jstep)); % Here we define where the top of the spatter pile is based on the number of clasts
    %xpiletopnode = 259; 

% System is separated into block matrices
% Inside each block matrix (ground, clast, and air) we calculate heat
% transfer by conduction only
% dT/dt=kappa*d^2T/dx^2
% Use forward difference in time for the lefthand side:
%   (Tnew(i) - Ti(i))/dt = dT/dt
% Use central differences in space for the righthand side:
%   (Tnew(i+1) -2*Tnew(i) + Tnew(i-1))/(dx^2) = d^2T/dx^2
% Combine to get:
%   (Tnew(i) - Ti(i))/dt = kappa*((Tnew(i+1) -2*Tnew(i) + Tnew(i-1))/(dx^2))
% Set equal to Ti(i)/dt to get:
% -Ti(i)/dt = (kappa*((Tnew(i+1) -2*Tnew(i) + Tnew(i-1))/(dx^2))) - Tnew(i)/dt
% Separate out all the variables
% Ti(i)/dt = -kappa*Tnew(i-1)/(dx^2) + kappa*2*Tnew(i)/(dx^2) + Tnew(i)/dt
% - kappa*Tnew(i+1)/(dx^2)
% Split it up how it will look in the matrix (remember that the LHS will be
% multplied by the Tnew at that location and is set equal to Ti(i) (RHS))
% L(i,i-1) = -kappaground/(dxstep^2)
% L(i,i) = 2*kappaground/(dxstep^2) + 1/dt
% L(i,i+1) = -kappaground/(dxstep^2)
% R(i,1) = Ti(i)/dt
% What we are doing here is using Ti to solve for our Tnew
% This gets us through one timestep. Then we set all the Tnew values that
% we just solved for to Ti and we solve everything again to get new Tnew
% values
% So let's actually do this

% First block matrix: Ground matrix
% Tells model to calculated heat transfer within the ground
for i=2:1:ypilebottomnode-1 % First block matrix ground from 1 node up from bottom of ground to one node below the boundary between the ground and the clast
L(i,i-1) = -kappaground/(dystep^2);
L(i,i) = 2*kappaground/(dystep^2) + 1/dt;
L(i,i+1) = -kappaground/(dystep^2);
R(i,1) = T0imp(i)/dt;
end

% Now we need to make sure that the heat flux is conserved across the
% boudnary between the ground and the clast
% Set the heat flux going into the boundary from the ground equal to the
% heat flux coming out of the boundary and into the clast
% (Ti(i+1) - Ti(i))/dx = (Ti(i) - Ti(i-1))/dx
% Set equal to Ti(i+1)
% Ti(i+1) = 2*Ti(i) - Ti(i-1)

% Let's actually put this into the matrix (remember that everything on
% the LHS of the matrix is multipled by the Ti at that location)
%L(xpilebottomnode,xpilebottomnode-1) = -1;
%L(xpilebottomnode,xpilebottomnode) = 2;
%R(xpilebottomnode,1) = Ti(xpilebottomnode+1);

% Heat flux continuity using energy balance 
L(ypilebottomnode,ypilebottomnode-1) = -kappaclast/(dystep^2);
L(ypilebottomnode,ypilebottomnode) = 2*kappaclast/(dystep^2) + 1/dt;
L(ypilebottomnode,ypilebottomnode+1) = -kappaclast/(dystep^2);
R(ypilebottomnode,1) = T0imp(i)/dt;


% Second block matrix: Clast Matrix CONDUCTION ONLY
% Tells model to calculate heat transfer within the clast
if jstep == 1
    for i=ypilebottomnode+1:1:yclastonetopnode-1 % Second block matrix from one node above the  boundary between the gound and the clast to one node below the boundary between the top of the spatter pile and the air
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end
    
 L(yclastonetopnode,yclastonetopnode-1) = (kappaclast)/(dystep^2);
    L(yclastonetopnode,yclastonetopnode) = ((-kappaclast)/(dystep^2)) + ((h)/(rhoclast*cpclast*dystep)) + (1/dt);
    R(yclastonetopnode,1) = T0imp(i)/dt + ((Tair*h)/(rhoclast*cpclast*dystep)) + ((-emis*sb*(T0imp(i)^4))/(rhoclast*cpclast*dystep)) + ((emis*sb*(Tair^4))/(rhoclast*cpclast*dystep)); 
    
% Third block matrix: Air Matrix
% Tells model to calculate heat transfer within the air
for i=yclastonetopnode+1:1:dysystemnode-1
L(i,i-1) = -kappaair/(dystep^2);
L(i,i) = 2*kappaair/(dystep^2) + 1/dt;
L(i,i+1) = -kappaair/(dystep^2);
R(i,1) = T0imp(i)/dt;
end    



% Define the last point in the matrix and the last point in the vector
L(dysystemnode,dysystemnode) = 1;
R(dysystemnode,1) = Tair;
end

if jstep == 2
   
    for i=ypilebottomnode+1:1:yclastonetopnode-1 % Second block matrix from one node above the  boundary between the gound and the clast to one node below the top of the first clast. If there is only one clast, this will be overwritten. But if there are more than one clast this is important for heat flux continuity across the clast/clast boundary
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end

    L(yclastonetopnode,yclastonetopnode-1) = -kappaclast/(dystep^2);
    L(yclastonetopnode,yclastonetopnode) = 2*kappaclast/(dystep^2) + 1/dt;
    L(yclastonetopnode,yclastonetopnode+1) = -kappaclast/(dystep^2);
    R(yclastonetopnode,1) = T0imp(i)/dt;

    for i=yclastonetopnode+1:1:yclasttwotopnode-1 % Second clast matrix from one node above the top of the first clast to one node below the top of the second clast
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end
    
 L(yclasttwotopnode,yclasttwotopnode-1) = (kappaclast)/(dystep^2);
    L(yclasttwotopnode,yclasttwotopnode) = ((-kappaclast)/(dystep^2)) + ((h)/(rhoclast*cpclast*dystep)) + (1/dt);
    R(yclasttwotopnode,1) = T0imp(i)/dt + ((Tair*h)/(rhoclast*cpclast*dystep)) + ((-emis*sb*(T0imp(i)^4))/(rhoclast*cpclast*dystep)) + ((emis*sb*(Tair^4))/(rhoclast*cpclast*dystep));

% Third block matrix: Air Matrix
% Tells model to calculate heat transfer within the air
for i=yclasttwotopnode+1:1:dysystemnode-1
L(i,i-1) = -kappaair/(dystep^2);
L(i,i) = 2*kappaair/(dystep^2) + 1/dt;
L(i,i+1) = -kappaair/(dystep^2);
R(i,1) = T0imp(i)/dt;
end    

% Define the last point in the matrix and the last point in the vector
L(dysystemnode,dysystemnode) = 1;
R(dysystemnode,1) = Tair;
end

if jstep == 3
    for i=ypilebottomnode+1:1:yclastonetopnode-1 % Second block matrix from one node above the  boundary between the gound and the clast to one node below the top of the first clast. If there is only one clast, this will be overwritten. But if there are more than one clast this is important for heat flux continuity across the clast/clast boundary
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end

    L(yclastonetopnode,yclastonetopnode-1) = -kappaclast/(dystep^2);
    L(yclastonetopnode,yclastonetopnode) = 2*kappaclast/(dystep^2) + 1/dt;
    L(yclastonetopnode,yclastonetopnode+1) = -kappaclast/(dystep^2);
    R(yclastonetopnode,1) = T0imp(i)/dt;

    for i=yclastonetopnode+1:1:yclasttwotopnode-1 % Second clast matrix from one node above the top of the first clast to one node below the top of the second clast
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end
    
    L(yclasttwotopnode,yclasttwotopnode-1) = -kappaclast/(dystep^2);
    L(yclasttwotopnode,yclasttwotopnode) = 2*kappaclast/(dystep^2) + 1/dt;
    L(yclasttwotopnode,yclasttwotopnode+1) = -kappaclast/(dystep^2);
    R(yclasttwotopnode,1) = T0imp(i)/dt;
    
    for i=yclasttwotopnode+1:1:yclastthreetopnode-1 % Second clast matrix from one node above the top of the first clast to one node below the top of the second clast
    L(i,i-1) = -kappaclast/(dystep^2);
    L(i,i) = 2*kappaclast/(dystep^2) + 1/dt;
    L(i,i+1) = -kappaclast/(dystep^2);
    R(i,1) = T0imp(i)/dt;
    end
    
    L(yclastthreetopnode,yclastthreetopnode-1) = (kappaclast)/(dystep^2);
    L(yclastthreetopnode,yclastthreetopnode) = ((-kappaclast)/(dystep^2)) + ((h)/(rhoclast*cpclast*dystep)) + (1/dt);
    R(yclastthreetopnode,1) = T0imp(i)/dt + ((Tair*h)/(rhoclast*cpclast*dystep)) + ((-emis*sb*(T0imp(i)^4))/(rhoclast*cpclast*dystep)) + ((emis*sb*(Tair^4))/(rhoclast*cpclast*dystep));

    
    
% Third block matrix: Air Matrix
% Tells model to calculate heat transfer within the air
for i=yclastthreetopnode+1:1:dysystemnode-1
L(i,i-1) = -kappaair/(dystep^2);
L(i,i) = 2*kappaair/(dystep^2) + 1/dt;
L(i,i+1) = -kappaair/(dystep^2);
R(i,1) = T0imp(i)/dt;
end    

% Define the last point in the matrix and the last point in the vector
L(dysystemnode,dysystemnode) = 1;
R(dysystemnode,1) = Tair;
end


% Now we need to make sure that the heat flux is conserved across the
% boundary between the top of the spatter pile and the air. But it is a
% little more complicated than the previous boundary because here we also
% have convection at the boundary. So we need to add convection into our
% energy balance.
% Energy balance for the boundary between the top of the spatter pile and
% the air:
% -kclast*((Tnew(i) - Tnew(i-1))/dx) - h*(Tnew(i) - Tnew(i+1)) =
% rho*cp*(dx/2)*((Tnew(i) - Ti(i))/dt)
% Set equal to Ti(i)/dt to get:
% Ti(i)/dt = +kappaclast*(2/dx)*((Tnew(i) - Tnew(i-1))/dx) + (2/rho*cp*dx)*h*(Tnew(i)
% - Tnew(i+1)) + Tnew(i)/dt
% Separate out all the variables:
% Ti(i)/dt = -kappaclast*(2/dx)*(Tnew(i-1)/dx) +
% kappaclast*(2/dx)*(Tnew(i)/dx) + (2/rho*cp*dx)*h*(Tnew(i)) + Tnew(i)/dt - (2/rho*cp*dx)*h*(Tnew(i+1))
% Let's actually put this into the matrix (remember that everything on
% the LHS of the matrix is multipled by the Ti at that location)

% Below is for radiation
   % L(xpiletopnode,xpiletopnode-1) = (-kappaclast)/(dxstep^2);
   % L(xpiletopnode,xpiletopnode) = ((kappaclast)/(dxstep^2)) + ((h)/(rhoclast*cpclast*dxstep)) + (1/dt);
   % R(xpiletopnode,1) = Ti(i)/dt + ((293*h)/(rhoclast*cpclast*dxstep)) + ((-emis*sb*(Ti(i)^4))/(rhoclast*cpclast*dxstep)) + ((emis*sb*(293^4))/(rhoclast*cpclast*dxstep));
%end

% Below is for convection
%L(xpiletopnode,xpiletopnode-1) = (-kappaclast)/(dxstep^2);
%L(xpiletopnode,xpiletopnode) = ((kappaclast)/(dxstep^2)) + ((h)/(rhoclast*cpclast*dxstep)) + (1/dt);
%R(xpiletopnode,1) = Ti(i)/dt + ((293*h)/(rhoclast*cpclast*dxstep));

% Below is for conduction
%L(xpiletopnode,xpiletopnode-1) = -kclast/(rhoclast*cpclast*(dxstep^2));
%L(xpiletopnode,xpiletopnode) = (kclast/(rhoclast*cpclast*(dxstep^2))) + (kair/(rhoclast*cpclast*(dxstep^2))) + 1/dt;
%L(xpiletopnode,xpiletopnode+1) = -kair/(rhoclast*cpclast*(dxstep^2));
%R(xpiletopnode,1) = Ti(i)/dt;



% Third block matrix: Air Matrix
% Tells model to calculate heat transfer within the air
%for i=xpiletopnode+1:1:dxsystemnode-1
%L(i,i-1) = -kappaair/(dxstep^2);
%L(i,i) = 2*kappaair/(dxstep^2) + 1/dt;
%L(i,i+1) = -kappaair/(dxstep^2);
%R(i,1) = Ti(i)/dt;
%end    

% Define the last point in the matrix and the last point in the vector
%L(dxsystemnode,dxsystemnode) = 1;
%R(dxsystemnode,1) = Tair;

% Now let's actually solve this!
% Obtaining solution for implicit temperature profile
% We input values for the LHS and the RHS of the matrix. LHSx=RHS. Now we
% need to solve for x, which is our Tnew.
T1imp = L\R;


if(timesum>(tclasttwo-dt) && timesum<(tclasttwo) && jstep<clastnum)
     T1imp(yclastonetopnode+1:yclasttwotopnode)= Tclast;
end

if(timesum>(tclasttwo) && timesum<(tclasttwo+dt) && jstep<clastnum)
    jstep=jstep+1;
end

if(timesum>(tclastthree-dt) && timesum<(tclastthree) && jstep<clastnum)
      T1imp(yclasttwotopnode+1:yclastthreetopnode)= Tclast;
end

if(timesum>(tclastthree) && timesum<(tclastthree+dt) && jstep<clastnum)
    jstep=jstep+1;
end
    


%if (timesum == tclastthree)
 %   Tnew(xclasttwotopnode+1:xclastthreetopnode)= Tclast;
  %  jstep == 3;
%end

%if (timesum > tnewclast && timesum < tnewclast + dt && jstep < clastnum)
%Tnew(xclastonetopnode+1:xclasttwotopnode)= Tclast; % added a new clast
%tnew = tnewclast + tbetween; %update time for next clast
%xpiletopnode = xpilebottomnode + (dxclastnode*jstep);
%jstep = jstep + 1;
%end

% Let's plot this implicit solution so we can see how the spatter pile
% cools with time
% Open figure
% Plotting implicit solution 
%subplot(2,1,1);
%plot(T0imp,y,'k');
%axis([0.9*Tground 1.1*Tclast 0 dysystem]);
%title(['Implicit: step=',num2str(t),' time,sec=',num2str(timesum)]);
%ylabel('x, m');
%xlabel('Temperature, K');
% Plot cooling rate
%subplot (2,1,2);
%plot (timesum, T0imp(yclastonetopnode-1),'o'); %This is the "location" of the thermocouple, we want it to be one node below the top of the first clast
timevec(b)=timesum;
tempvec(b)=T0imp(yclastonetopnode-1);
b=b+1;
%hold on;
disp(dt);
disp(timesum);
%disp(jstep);
%disp(Ti(322));
%disp(Ti(321));
%disp(Ti(320));
%disp(Ti(301));
%disp(Ti(282));
%disp(Ti(281));
%disp(Ti(280));
%disp(Ti(261));
%disp(Ti(242));
%disp(Ti(241));
%disp(Ti(240));
%disp(Ti(221));
%disp(Ti(202));
%disp(Ti(201));
disp(T0imp(yclastonetopnode-1));
disp(T0imp(yclastonetopnode));
disp(T0imp(yclastonetopnode+1));
%disp(T0imp(yclasttwotopnode+1));

%axis([0 500 300 1410]);
%xlabel('Time, sec');
%ylabel('Temperature, K');
%text (10,10,'1025 sec - 2nd flow');

% Stop for 0.1 second
%drawnow % draw above figure without delay
%pause(0.1);

% Stop for 0.1 second
%drawnow % draw above figure without delay
%pause(0.1);

% Add time counter
timesum = timesum+dt; 

% reassign temperature profiles for next step. This is where we say, okay,
% we solved everything to get out new temperature (Tnew) using Ti. Now we
% say, let's make all the Tnew that we solved for be Ti, then we can solve
% everything again and get new Tnew. 
T0imp = T1imp;
%if Ti(dxclastnode-1) < 890.
 %   break
%end
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

%A = jstep;
%plot(A)
