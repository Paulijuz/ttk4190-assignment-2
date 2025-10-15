% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

% Add folder for 3-D visualization files
addpath(genpath('flypath3d_v2'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
T_final = 1000;	        % Final simulation time (s)
h = 0.1;                % Sampling time (s)

psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_ref   = 7;            % desired surge speed (m/s)

% initial states
eta_0 = [0 0 0]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
x = [nu_0' eta_0' delta_0 n_0]'; % The state vector can be extended with addional states here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:h:T_final;                % Time vector
nTimeSteps = length(t);         % Number of time steps

simdata = zeros(nTimeSteps, 13); % Pre-allocate matrix for efficiency

for i = 1:nTimeSteps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1a) Add current here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uc = 0;
    vc = 0;
    nu_c = [ uc vc 0 ]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1c) Add wind here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ywind = 0;
    Nwind = 0;
    tau_wind = [0 Ywind Nwind]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add a reference model here 
    % Define it as a function
    % check eq. (15.143) in (Fossen, 2021) for help
    %
    % The result should look like this:
    % xd_dot = ref_model(xd,psi_ref(i));
    % psi_d = xd(1);
    % r_d = xd(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi_d = psi_ref;
    r_d = 0;
    u_d = U_ref;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add the heading controller here 
    % Define it as a function
    %
    % The result should look like this:
    % delta_c = PID_heading(e_psi,e_r,e_int);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta_c = 0.1;              % rudder angle command (rad)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3e) Add open loop speed control here
    % Define it as a function
    %
    % The result should look like this:
    % n_c = open_loop_speed_control(U_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_c = 10;                   % propeller speed [radians per second (rps)]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3f) Replace the open loop speed controller, 
    % with a closed loop speed controller here 
    % Define it as a function
    %
    % The result should look like this:
    % n_c = closed_loop_speed_control(u_d,e_u,e_int_u);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,tau_total] = ship(x,u,nu_c,tau_wind);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d];     
 
    % Euler integration
    % x = euler2(xdot,x,h); 
    % Runge Kutta 4 integration
    x = rk4(@ship,h,x,u,nu_c,tau_wind);

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end

end
simdata = simdata(1:i,:);
t = t(1:i);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
x           = simdata(:,4);                 % m
y           = simdata(:,5);                 % m
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c         = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         =  simdata(:,13);               % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
%%
figure(3)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi_deg,t,psi_d_deg,'linewidth',2);
title('Actual and desired yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 
legend('actual yaw','desired yaw')
subplot(313)
plot(t,r_deg,t,r_d_deg,'linewidth',2);
title('Actual and desired yaw rates'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 
legend('actual yaw rate','desired yaw rate')

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')
%% Create objects for 3-D visualization 
% Since we only simulate 3-DOF we need to construct zero arrays for the 
% excluded dimensions, including height, roll and pitch
z = zeros(length(x),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);

% create object 1: ship (ship1.mat)
new_object('flypath3d_v2/ship1.mat',[x,y,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(x)),max(abs(y)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);

% Plot trajectories 
flypath('flypath3d_v2/ship1.mat',...
'animate','on','step',500,...
'axis','on','axiscolor',[0 0 0],'color',[1 1 1],...
'font','Georgia','fontsize',12,...
'view',[-25 35],'window',[900 900],...
'xlim', [min(y)-0.1*max(abs(y)),max(y)+0.1*max(abs(y))],... 
'ylim', [min(x)-0.1*max(abs(x)),max(x)+0.1*max(abs(x))], ...
'zlim', [-max(max(abs(x)),max(abs(y)))/100,max(max(abs(x)),max(abs(y)))/20]); 