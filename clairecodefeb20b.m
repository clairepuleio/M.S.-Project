% Solution of 2D temperature equation on a regular grid 
% for a non-moving medium with constant conductivity;
% comparison of implicit and explicit method


% Clean all variables
clear all;
% Clear all figures
clf;

% Numerical model parameters
% Boundary conditions = constant temperature
% setup correspond to a squre temperature wave
% Model size, m
dxsystem=10; % Horizontal
dysystem=10; % Vertical
yclastonetop=5;

% Numbers of nodes
dxsystemnode=11; % Horizontal
dysystemnode=11; % Vertical
yclastonetopnode=6;
% Grid step
xstp=dxsystem/(dxsystemnode-1); % Horizontal
ystp=dysystem/(dysystemnode-1); % Vertical
tnum=10; % number of timesteps

% Material properties
kclast=1.5; % Thermal conductivity, W/m/K
cpclast=1490; % Heat capacity, J/kg/K
rhoclast=1300; % Density, kg/m^3
kappaclast=kclast/rhoclast/cpclast; % Thermal diffusivity, m^2/s

kair=0.0264;
cpair=1006.;
rhoair=1.204; 
kappaair=kair/rhoair/cpair; 

% Making vectors for nodal points positions (basic nodes)
x=0:xstp:dxsystem; % Horizontal
y=0:ystp:dysystem; % Vertical

% Timestep
dtexp=min(xstp,ystp)^2/3/kappaclast; % Limitation for explicit timestep 
dt=dtexp/100000; % Timestep, s


% Define initial temperature profile
Tair=300; % background temperature, K
Tclast=1400; % temperature wave, K
% Creatiing arrays
%t0exp=tback*ones(ynum,xnum); % for explicit solving 
 % for implicit solving
t0imp=Tair*ones(dysystemnode,dxsystemnode); % for implicit solving
for i=1:1:dxsystemnode
    for j=1:1:dysystemnode
        % Temperature wave
        if(i>1 && i<dxsystemnode && j>1 && j<=yclastonetopnode) 
           % t0exp(i,j)=twave; % for explicit solving
            t0imp(i,j)=Tclast; % for implicit solving
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
    for i=1:1:dxsystemnode
        for j=1:1:dysystemnode
            % Global index
            k=(j-1)*dxsystemnode+i;
            % Boundary nodes
           
   % Bottom Boundary
                if(j==1)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tclast; 
                   
                end
                % Top boundary
                if(j==dysystemnode)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tair;   
                    
                end
                % Left boundary lower
                if(i==1 && j>1 && j<=yclastonetopnode)
                    % Constant temperature: T(i,j)=tback
                    L(k,k)=1;
                    R(k,1)=Tclast;     
                    
                end
                % Left boundary higher
                if(i==1 && j>yclastonetopnode && j<dysystemnode)
                    L(k,k)=1;
                    R(k,1)=Tair;
                    
                end
                % Right boundary lower
                if(i==dxsystemnode && j>1 && j<=yclastonetopnode)
                   % Constant temperature: T(i,j)=Tclast
                    L(k,k)=1;
                    R(k,1)=Tclast;      
                    
                end
                % Right boundary higher
                if(i==dxsystemnode && j>yclastonetopnode && j<dysystemnode)
                    L(k,k)=1;
                    R(k,1)=Tair;
                    
                end 
               % Clast top boundary
                 if(i>1 && i<dxsystemnode && j==yclastonetopnode)
                    L(k,k-dxsystemnode)=kclast/((ystp^2)*rhoclast*cpclast); % coefficient for T(i,j-1) OR (k,k-N) where N=dxsystemnode
                    L(k,k+dxsystemnode)=kclast/((ystp^2)*rhoclast*cpclast); % coefficient for T(i,j+1) OR (k,k+N) where N=dxsystemnode
               %     L(k,k-1)=kclast/((xstp^2)*rhoclast*cpclast);
               %     L(k,k+1)=kclast/((xstp^2)*rhoclast*cpclast);
                    L(k,k)=((-kclast-kclast)/((ystp^2)*rhoclast*cpclast))+(1/dt); % coefficient for T(i,j)
                    R(k,1)=t0imp(i,j)/dt;
                    
                    
                 end
                 
        
            % Internal nodes clast
          
        
          
         if(i>1 && i<dxsystemnode && j>1 && j<yclastonetopnode)
                L(k,k-dxsystemnode)=-kappaclast/ystp^2; % coefficient for T(i,j-1)
                L(k,k+dxsystemnode)=-kappaclast/ystp^2; % coefficient for T(i,j+1)
                L(k,k-1)=-kappaclast/xstp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-kappaclast/xstp^2; % coefficient for T(i+1,j)
                L(k,k)=1/dt+2*kappaclast/ystp^2+2*kappaclast/xstp^2; % coefficient for T(i,j)
                R(k,1)=t0imp(i,j)/dt;
            %   disp('calculating clast nodes')
         end
            % Internal nodes air
         if(i>1 && i<dxsystemnode && j>yclastonetopnode && j<dysystemnode)
                L(k,k-dxsystemnode)=-kappaclast/ystp^2; % coefficient for T(i,j-1)
                L(k,k+dxsystemnode)=-kappaclast/ystp^2; % coefficient for T(i,j+1)
                L(k,k-1)=-kappaclast/xstp^2; % coefficient for T(i-1,j)
                L(k,k+1)=-kappaclast/xstp^2; % coefficient for T(i+1,j)
                L(k,k)=1/dt+2*kappaclast/ystp^2+2*kappaclast/xstp^2; % coefficient for T(i,j)
                R(k,1)=t0imp(i,j)/dt;
            %    disp('calc air nodes')
         end
        end
    end
    % Obtaining solution
    S=L\R;
    % Reloading solution to the new temperature array
    for i=1:1:dxsystemnode
        for j=1:1:dysystemnode
            % Global index
            k=(j-1)*dxsystemnode+i;
            % Reload solution
            t1imp(i,j)=S(k);
        end
    end
    


%  spy(L)

display(t0imp(50))
display(t0imp(61))
display(t0imp(72))


      
    % Open figure
    figure(1);
    % Plotting implicit solution
    subplot(1,1,1);
    pcolor(x,y,t1imp);
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
    t0imp=t1imp;
    end
        
