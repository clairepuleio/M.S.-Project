% 2D 1 clast conduction and convection
% Clear all variables
clear all;
% Clear all figures
clf;
 
%CHANGE THIS ABOUT CODE EACH TIME
dxclast=0.05; % width of clast (m)
tnum=500000; % number of timesteps
clastnum=1; % number of clasts in spatter pile
%tclasttwo=5; % time clast two is deposited
%tclastthree=10; % time clast three is deposited
 
% xtracker is where horizontally in our clast we want to track the temp

xtracker=151;
 
% Numerical model parameters
 
dxsystem=3; % Horizontal model size (m)
dysystem=3; % Vertical model size (m)
dxsystemnode=301; % Horizontal model size (nodes)
dysystemnode=301; % Vertical model size (nodes)
xstp=dxsystem/(dxsystemnode-1); % Horizontal grid step
ystp=dysystem/(dysystemnode-1); % Vertical grid step
 
dxclastnode=dxclast/xstp; % clast size (nodes)
 
ypilebottom=1.; % location of bottom of spatter pile (m)
yclastonetop=ypilebottom+dxclast; % location of top of clast one (m)
yclasttwotop=yclastonetop+dxclast; % location of top of clast two (m)
yclastthreetop=yclasttwotop+dxclast; % location of top of clast three (m)
xclastleft=1.48; % location of left side of spatter pile (m)
xclastright=xclastleft+dxclast; % location of right side of spatter pile (m)
 
ypilebottomnode=101; % location of bottom of spatter pile (nodes)
yclastonetopnode=ypilebottomnode+dxclastnode; % location of top of clast one (nodes)
yclasttwotopnode=yclastonetopnode+dxclastnode; % location of top of clast two (ndoes)
yclastthreetopnode=yclasttwotopnode+dxclastnode; % location of top of clast three (nodes)
xclastleftnode=149; % location of left side of spatter pile (nodes)
xclastrightnode=(xclastleftnode+dxclastnode)-1; % location of right side of spatter pile (nodes)
 
dxpile=dxclast*clastnum; % thickness of spatter pile
 
jstep=1; %step after every new clast added to spatter pile
 
% Material properties
kclast=0.3601; % Thermal conductivity clast/ground, W/m/K
cpclast=1490.; % Heat capacity clast/ground, J/kg/K
rhoclast=850.; % Density clast/ground, kg/m^3
kappaclast=kclast/rhoclast/cpclast; % Thermal diffusivity clast/ground, m^2/s
 
kair=0.0264; % Thermal conductivity air, W/m/K
cpair=1006.; % Heat capacity air, J/kg/K
rhoair=1.204; % Density air, kg/m^3
kappaair=kair/rhoair/cpair; % thermal diffusivity air, m^2/s
 
h=8.; % Heat transfer coefficient, W/m^2 K
emis = 0.96; %emissivity taken from engineering toolbox SUBJECT TO CHANGE
sb = 0.0000000567; %stefan-boltzmann constant
 
% Making vectors for nodal points positions (basic nodes)
x=0:xstp:dxsystem; % Horizontal
y=0:ystp:dysystem; % Vertical
 
% Timestep
dtexp=min(xstp,ystp)^2/3/kappaclast; % Limitation for explicit timestep 
dt=dtexp/10; % Timestep, s
  
% Define initial temperature profile
Tair=293.; % background temperature, K
Tclast=1373.; % temperature wave, K
Tground=293.;
% Creatiing arrays
%t0exp=tback*ones(ynum,xnum); % for explicit solving 
 % for implicit solving
T0imp=Tair*ones(dysystemnode,dxsystemnode); % for implicit solving
for i=1:1:dysystemnode
    for j=1:1:dxsystemnode
        % Temperature wave
        if((i)>=ypilebottomnode+1 && (i)<=yclastonetopnode && (j)>=xclastleftnode && (j)<=xclastrightnode) 
           % t0exp(i,j)=twave; % for explicit solving
            T0imp(i,j)=Tclast; % for implicit solving
        end
    end
end
 
% Time cycle
timesum=0; % Elapsed time
for t=1:1:tnum
 
    % Matrix of coefficients initialization for implicit solving
    L=sparse(dxsystemnode*dysystemnode,dxsystemnode*dysystemnode);
    % Vector of right part initialization for implicit solving
    R=zeros(dxsystemnode*dysystemnode,1);
    
    % Implicit solving of 2D temperature equation:
    % dT/dt=kappa*(d2T/dx2+d2T/dy2)
    % Composing matrix of coefficients L()
    % and vector (column) of right parts R()
    % Process all grid points
    for i=1:1:dysystemnode
        for j=1:1:dxsystemnode
            % Global index
            k=(j-1)*dysystemnode+i;
            % Boundary nodes
           
   % Bottom Boundary
                if(i==1)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tground; 
                   
                end
                % Top boundary
                if(i==dysystemnode)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tair;   
                    
                end
                % Left boundary ground
                if(j==1 && i>1 && i<=ypilebottomnode)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tground;       
                end
               
        % Left boundary air
        if(j==1 && i>ypilebottomnode && i<dysystemnode)
            L(k,k)=1;
            R(k,1)=Tair;
        end
    
                % Right boundary ground
                if(j==dxsystemnode && i>1 && i<=ypilebottomnode)
                   % Constant temperature: T(i,j)=Tclast
                    L(k,k)=1;
                    R(k,1)=Tground;      
                    
                end
               
        % Right boundary air
        if(j==dxsystemnode && i>ypilebottomnode && i<dysystemnode)
            L(k,k)=1;
            R(k,1)=Tair;
        end
        
% One clast in pile
        if jstep == 1
 
           % Clast bottom boundary
         if(j>xclastleftnode && j<xclastrightnode && i==ypilebottomnode+1)
            L(k,k-1)=-kclast/(rhoclast*cpclast*(ystp^2));
            L(k,k)=((kclast+kclast)/(rhoclast*cpclast*(ystp^2)))+(1/dt);
            L(k,k+1)=-kclast/(rhoclast*cpclast*(ystp^2));
            R(k,1)=T0imp(i,j)/dt;
        end
 
               % Clast top boundary
        if(j>xclastleftnode && j<xclastrightnode && i==yclastonetopnode)
            L(k,k-1)=-kclast/((ystp^2)*rhoclast*cpclast);
            L(k,k)=((kclast)/((ystp^2)*rhoclast*cpclast))+(h/(rhoclast*cpclast*ystp))+(1/dt);
            R(k,1)=(T0imp(i,j)/dt)+((Tair*h)/(rhoclast*cpclast*ystp));
        end
                 
               % Clast left boundary 
        if(j==xclastleftnode && i>=ypilebottomnode+1 && i<=yclastonetopnode)   
            L(k,k)=((kclast)/((xstp^2)*rhoclast*cpclast))+(h/(rhoclast*cpclast*xstp))+(1/dt);
            L(k,k+dysystemnode)=((-kclast)/((xstp^2)*rhoclast*cpclast));
            R(k,1)=(T0imp(i,j)/dt)+((Tair*h)/(rhoclast*cpclast*xstp));
        end   
          
               % Clast right boundary
        if(j==xclastrightnode && i>=ypilebottomnode+1 && i<=yclastonetopnode)
            L(k,k-dysystemnode)=((-kclast)/((xstp^2)*rhoclast*cpclast)); 
            L(k,k)=((kclast)/((xstp^2)*rhoclast*cpclast))+(h/(rhoclast*cpclast*xstp))+(1/dt); 
           R(k,1)=(T0imp(i,j)/dt)+((Tair*h)/(rhoclast*cpclast*xstp));
        end 
        
               % Left ground/air boundary
                if(j>1 && j<xclastleftnode && i==ypilebottomnode+1)   
            L(k,k-1)=-kclast/(rhoclast*cpclast*(ystp^2));
            L(k,k)=((kclast+kair)/(rhoclast*cpclast*(ystp^2)))+(1/dt);
            L(k,k+1)=-kair/(rhoclast*cpclast*(ystp^2));
            R(k,1)=T0imp(i,j)/dt;
               end
               
               % Right ground/air boundary
                if(j>xclastrightnode && j<dxsystemnode && i==ypilebottomnode+1)
                  L(k,k-1)=-kclast/(rhoclast*cpclast*(ystp^2));
                  L(k,k)=((kclast+kair)/(rhoclast*cpclast*(ystp^2)))+(1/dt);
                  L(k,k+1)=-kair/(rhoclast*cpclast*(ystp^2));
                  R(k,1)=T0imp(i,j)/dt;
               end
                
        % Internal nodes ground
     if(j>1 && j<dxsystemnode && i>1 && i<=ypilebottomnode)
        L(k,k-dysystemnode)=-kappaclast/xstp^2;
        L(k,k+dysystemnode)=-kappaclast/xstp^2;
        L(k,k-1)=-kappaclast/ystp^2;
        L(k,k+1)=-kappaclast/ystp^2;
        L(k,k)=1/dt+2*kappaclast/ystp^2+2*kappaclast/xstp^2;
        R(k,1)=T0imp(i,j)/dt;
    end        
            % Internal nodes clast
          
         if(j>xclastleftnode && j<xclastrightnode && i>ypilebottomnode+1 && i<yclastonetopnode)
                L(k,k-dysystemnode)=-kappaclast/xstp^2; % coefficient for T(i,j-1)
                L(k,k+dysystemnode)=-kappaclast/xstp^2; % coefficient for T(i,j+1)
                L(k,k-1)=-kappaclast/ystp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-kappaclast/ystp^2; % coefficient for T(i+1,j)
                L(k,k)=1/dt+2*kappaclast/ystp^2+2*kappaclast/xstp^2; % coefficient for T(i,j)
                R(k,1)=T0imp(i,j)/dt;
            %   disp('calculating clast nodes')
         end
         
            % Internal nodes left air
            if(j>1 && j<xclastleftnode && i>ypilebottomnode && i<=yclastonetopnode)
               L(k,k-dysystemnode)=-kappaair/(xstp^2); % coefficient for T(i,j-1)
                L(k,k+dysystemnode)=-kappaair/(xstp^2); % coefficient for T(i,j+1)
                L(k,k-1)=-kappaair/(ystp^2); % coefficient for T(i-1,j)
                L(k,k+1)=-kappaair/(ystp^2); % coefficient for T(i+1,j)
                L(k,k)=(1/dt)+(2*kappaair/(ystp^2))+(2*kappaair/(xstp^2)); % coefficient for T(i,j)
                R(k,1)=T0imp(i,j)/dt;
            end
            % Internal nodes right air
             if(j>xclastrightnode && j<dxsystemnode && i>ypilebottomnode && i<=yclastonetopnode)
               L(k,k-dysystemnode)=-kappaair/(xstp^2); % coefficient for T(i,j-1)
                L(k,k+dysystemnode)=-kappaair/(xstp^2); % coefficient for T(i,j+1)
                L(k,k-1)=-kappaair/(ystp^2); % coefficient for T(i-1,j)
                L(k,k+1)=-kappaair/(ystp^2); % coefficient for T(i+1,j)
                L(k,k)=(1/dt)+(2*kappaair/(ystp^2))+(2*kappaair/(xstp^2)); % coefficient for T(i,j)
                R(k,1)=T0imp(i,j)/dt;
            end
            % Internal nodes upper air
         if(j>1 && j<dxsystemnode && i>yclastonetopnode && i<dysystemnode)
                L(k,k-dysystemnode)=-kappaair/xstp^2; % coefficient for T(i,j-1)
                L(k,k+dysystemnode)=-kappaair/xstp^2; % coefficient for T(i,j+1)
                L(k,k-1)=-kappaair/ystp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-kappaair/ystp^2; % coefficient for T(i+1,j)
                L(k,k)=1/dt+2*kappaair/ystp^2+2*kappaair/xstp^2; % coefficient for T(i,j)
                R(k,1)=T0imp(i,j)/dt;
            %    disp('calc air nodes')
         end
        end
  
        end
    end
    % Obtaining solution
    S=L\R;
    % Reloading solution to the new temperature array
    for i=1:1:dysystemnode
        for j=1:1:dxsystemnode
            % Global index
            k=(j-1)*dysystemnode+i;
            % Reload solution
            T1imp(i,j)=S(k);
        end
    end
 
% Adding in clasts
%if(timesum>(tclasttwo-dt) && timesum<(tclasttwo) && jstep<clastnum)
%     T1imp(yclastonetopnode:yclasttwotopnode,xclastleftnode:xclastrightnode)=Tclast;
%end
%if(timesum>(tclasttwo) && timesum<(tclasttwo+dt) && jstep<clastnum)
%    jstep=jstep+1;
%end
 
%if(timesum>(tclastthree-dt) && timesum<(tclastthree) && jstep<clastnum)
%       T1imp(yclasttwotopnode:yclastthreetopnode,xclastleftnode:xclastrightnode)=Tclast;
%end
 
%if(timesum>(tclastthree) && timesum<(tclastthree+dt) && jstep<clastnum)
%    jstep=jstep+1;
%end
    
 display(timesum)
display(T0imp(yclastonetopnode-1,xtracker))
     % Open figure
    figure(1);
    % Plotting implicit solution
    subplot(1,1,1);
    pcolor(x,y,T0imp);
    shading interp;
    colorbar;
    axis ij image;
  %     caxis([tback-100 twave+100]);
    title(['Implicit: step=',num2str(t),' time,Myr=',num2str(timesum/(1e+6*365.25*24*3600))]);
   
    % Stop for .1 second
    pause(.1);
    
    % Add time counter
    timesum=timesum+dt;
    % Reassign temperature profiles for the next step
   % t0exp=t1exp;
 
if T0imp(yclastonetopnode-1,xtracker) < 895.5
break
end
 
   
T0imp=T1imp;
end
