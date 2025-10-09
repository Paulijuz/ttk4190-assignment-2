function xdot = ship(x,u)
% xdot = ship(x,u) returns the time derivative of the state vector: 
% x = [ u v r x y psi delta n ]' for a ship with L = 161 m where:
%
% u     = surge velocity, must be positive  (m/s)    
% v     = sway velocity                     (m/s)
% r     = yaw velocity                      (rad/s)
% x     = position in x-direction           (m)
% y     = position in y-direction           (m)
% psi   = yaw angle                         (rad)
% delta = actual rudder angle               (rad)
% n     = actual shaft velocity             (rps)
% 
% The input vector is :
%
% u       = [ delta_c  n_c ]'  where
%
% delta_c = commanded rudder angle          (rad)
% n_c     = commanded shaft velocity        (rpm)
%
% Author:    name
% Date:      date

% Check of input and state dimensions
if (length(x)~= 8),error('x-vector must have dimension 8 !');end
if (length(u)~= 2),error('u-vector must have dimension 2 !');end

% Dimensional states and input
delta_c = u(1); 
n_c     = u(2);

nu    = x(1:3);
eta   = x(4:6);
delta = x(7);
n     = x(8); 

% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia (kg m^2)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
KT = 0.7;               % propeller coefficient (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (m/s^3)

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add added mass here

MA = -[
    Xudot,  0,      0;
    0,      Yvdot,  Yrdot;
    0,      Nvdot,  Nrdot;
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];

Minv = invQR(MRB + MA);

% input matrix
t_thr = 0.05;           % thrust deduction number
X_delta2 = 0;         % rudder coefficients (Section 9.5)
Y_delta = 0;      
N_delta = -1e7;
Bi = [ (1-t_thr) X_delta2
        0        Y_delta
        0        N_delta  ];
    
% state-dependent time-varying matrices
CRB = m * nu(3) * [ 0 -1 -xg 
                    1  0  0 
                    xg 0  0  ];
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Coriolis due to added mass here

CA = m2c(MA, nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add linear damping here

T1 = 20;
T2 = 20;
T6 = 10;

D11 = (m - Xudot) / T1;
D22 = (m - Yvdot) / T2;
D33 = (Iz - Nrdot) / T6;

D = diag([D11, D22, D33]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add nonlinear damping here

Cd_2D = Hoerner(B,T);

u = x(1);
v = x(2);
r = x(3);

Ucf = 0;
Ycf = 0;
Ncf = 0;

% Strip theory: cross−flow drag integrals
dx = L/10; % 10 strips
for xL = -L/2:dx:L/2
    Ucf = abs(v + xL * r) * (v + xL * r);
    Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx; % sway force
    Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx; % yaw moment
end

Rn = L/(10^-6)*abs(u);

epsilon = 0.001;

Cf = 0.075 / ((log10(Rn)-2)^2 + epsilon);

S = 2*T*B + 2*T*L + L*B;

k = 0.1;

Xcf = - 0.5 * rho * S * (1 + k) * Cf * abs(u) * u;

d = -[
    Xcf;
    Ycf;
    Ncf;
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = Rzyx(0,0,eta(3));

% thrust 
thr = rho * Dia^4 * KT * abs(n) * n;    % thrust command (N)

% ship dynamics
u = [ thr delta ]';
tau = Bi * u;

% Denne linjen har blitt endret på fra den utdelte koden.
% La til "+ Ca + D" og "- d".
nu_dot = Minv * (tau - (CRB + CA + D) * nu - d);

eta_dot = R * nu;    

% Rudder saturation and dynamics (Sections 9.5.2)
if abs(delta_c) >= delta_max
    delta_c = sign(delta_c)*delta_max;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max
    delta_dot = sign(delta_dot)*Ddelta_max;
end    

% propeller dynamics
n_dot = (1/10) * (n_c - n);

xdot = [nu_dot' eta_dot' delta_dot n_dot]';

end