%% Deployment Scheme: Circular LCV AERO 449 Martin Kamme

% Circular LCV Orbit Below Operational Orbit

%% 
close all
clear all
clc

mumoon = 4903;                              % Moon's Gravitational Parameter km^3/s^2
T = 14/(1000*1000);                         % Thrust [kg*km/s^2]
Isp = 1100;                                 % Isp [seconds]
m = 179;
m0 = 179;                                   % Mass [kg]
g0 = .00162;                                % km/s^2

% Operational Orbit--------------------------------------------------------
inc_op = deg2rad(55);                       % Inclination [rad]
ecc_op = 0;                                 % Eccentricity
RAAN_op = 0;                                % Right Ascension of Ascending Node [rad]
a_op = 7500;                                % Semi Major Axis [km]   
T_op = ((2*pi)/sqrt(mumoon))*(a_op^(3/2));  % Period [seconds]
V_op = sqrt(mumoon/a_op);                   % Velocity [km/s]
h_op = sqrt(a_op*mumoon);                   % Angular momentum [km^2/s]
theta_op = 0;                               % True anomaly [rad]
per_op = 0;                                 % Argument of perilune [rad]
n_op = (2*pi)/T_op;                         % Mean Motion [rad/s]

[R_op,V_op] = coes2rvKAMME(h_op,inc_op,RAAN_op,ecc_op,per_op,theta_op,mumoon);



% LCV phaser Orbit---------------------------------------------------------
below_opTEST = linspace(75,170,20);

 for ii = 1:length(below_opTEST)
below_op = below_opTEST(ii);                % LCV alt below Operational orbit [km]
%below_op = 120;

inc_lcv = deg2rad(55);                      % Inclination [inc]
ecc_lcv = 0;                                % Eccentricity
RAAN_lcv = 0;                               % Right Ascension of Ascending Node [rad]
r_lcv = a_op - below_op;                    % Orbital Radius [km]
T_lcv = ((2*pi)/sqrt(mumoon))*(r_lcv^(3/2));% Period [seconds]
V_lcv = sqrt(mumoon/r_lcv);                   % Velocity [km/s]
h_lcv = sqrt(r_lcv*mumoon*(1-ecc_lcv^2));   % Angular Momentum [km^2/s]
theta_lcv = deg2rad(0);                     % True anomaly [rad]
per_lcv = deg2rad(0);                       % Argument of perilune [rad] define as being located 0 degrees from the ascending node
n_lcv = 2*pi/(T_lcv);                       % Mean Motion [rad/s] 
thetaINT = deg2rad(270);

LCVplotRadius(ii) = r_lcv;

[R_lcv,V_lcv] = coes2rvKAMME(h_lcv,inc_lcv,RAAN_lcv,ecc_lcv,per_lcv,theta_lcv,mumoon);


% State Vectors and Things-------------------------------------------------
[R_lss,V_lss] = coes2rvKAMME(h_op,inc_op,RAAN_op,ecc_op,per_op,thetaINT,mumoon);
[R_int,V_int] = coes2rvKAMME(h_op,inc_op,RAAN_op,ecc_op,per_op,thetaINT,mumoon);

state_op = [R_op;V_op];
state_lcv = [R_lcv;V_lcv];
state_lss = [R_lss;V_lss];
state_NonImpulsive = [R_lcv;V_lcv;m];

time = T_op;

time_lcv = T_lcv;
tspan = [0 time];
tspan_lcv = [0 time_lcv];


options = odeset('RelTol',1e-8,'AbsTol', 1e-8);
[tnew_op,orbit_op] = ode45(@orbitPROP,tspan,state_op,options,mumoon);
[tnew_lcv,orbit_lcv] = ode45(@orbitPROP,tspan_lcv,state_lcv,options,mumoon);
%[tnew_lss,orbit_lss] = ode45(@orbitPROP,tspan_lss,state_lss,options,mumoon);


% Attempt the non impulsive with Encke's
dr = [0;0;0];
dv = [0;0;0];
dt = 1; % set time step as 1 second

rpert = R_lcv;
vpert = V_lcv;
rosc = R_lcv;
vosc = V_lcv;

radius = 12;
t = 0; % initiliaze time
i = 1; % counter for tracking values
m = 179; 

while radius <= 7500
    [rosc,vosc] = keplerUniversal(rosc,vosc,dt,mumoon);
    dmdt = -T/(Isp*g0);
    m = dmdt*dt + m;
    state = [rosc;vosc;m];
    ap = nonImpulsive(state);
    da = ap;
    dv = da*dt + dv;
    dr = (1/2)*da*dt^2 + dv*dt + dr;
    rpert = rosc + dr;
    vpert = vosc + dv;
    t = t + dt;
    
    % store info for plotting
    rplot(:,i) = rpert.*1;
    vplot(:,i) = vpert.*1;
    radius = norm(rplot(:,i));
    
    % rectify every step
    rosc = rpert;
    vosc = vpert; 
    dr = [0;0;0];
    dv = [0;0;0];
    
    i = i+1;
    
end

TIME(ii) = t;
deltaV(ii) = 1000*abs(Isp*g0*log(1-(T/(m0*g0*Isp))*t));
time_lss = t;
tspan_NonImpulsive = [0 time_lss];
[t_NonImpulsive,orbit_NonImpulsive] = ode45(@orbitPropNonImpulsive,tspan_NonImpulsive,state_NonImpulsive,options,mumoon);


% Plot the things----------------------------------------------------------

[x,y,z] = sphere;
x = x*1731;
y = y*1731;
z = z*1731;

figure
hold on 
surf(x,y,z,'FaceColor',[0.8500 0.3250 0.0980])
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
plot3(orbit_op(:,1),orbit_op(:,2),orbit_op(:,3),'g')
plot3(orbit_lcv(:,1),orbit_lcv(:,2),orbit_lcv(:,3),'m')
%plot3(R_int(1),R_int(2),R_int(3),'r*')
%plot3(orbit_NonImpulsive(:,1),orbit_NonImpulsive(:,2),orbit_NonImpulsive(:,3),'r')
plot3(rplot(1,:),rplot(2,:),rplot(3,:),'b')
%plot3(rplot(1,end),rplot(2,end),rplot(3,end),'r*')
%plot3(R_lcv(1),R_lcv(2),R_lcv(3),'b*')

timeToWait(ii) = ((pi/4))/(n_lcv - n_op);

missionTimeLCV(ii) = (7*timeToWait(ii))/(3600*24);
missionTimeTotal(ii) = (7*timeToWait(ii) + TIME(ii))/(3600*24);
 end

TIMEplot = TIME./(3600);
massProp = TIME.*(T/(Isp*g0));

%% Plots
clc


figure
hold on
grid on
plot(below_opTEST,missionTimeLCV)
plot(below_opTEST,missionTimeTotal)
xlabel('Beginning Elevation Below Operational orbit')
ylabel('Time (days)')
legend('LCV Mission Time','Total Mission Time')


figure
hold on 
grid on
plot(below_opTEST,TIMEplot)
xlabel('Beginning Elevation Below Operational orbit')
ylabel('Time (hours)')
title('Thrust Time')

figure 
hold on 
grid on 
title('Propellant Mass for Deployment')
plot(below_opTEST,massProp)
xlabel('Beginning Elevation Below Operational orbit')
ylabel('Propellant Mass (kg)')

figure
hold on
grid on
title('DeltaV')
plot(below_opTEST,deltaV)
xlabel('Beginning Elevation Below Operational orbit')
ylabel('DeltaV (m/s)')
%% Functions

% Standard Functions
function [r,v] = keplerUniversal(r0,v0,t,mu)
%input vectors as columns

%Purpose:
%Most effecient way to propagate any type of two body orbit using kepler's
%equations.
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%t                      [1 x N]                         time vector in
%                                                       seconds
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
% Outputs:
%---------                                                                %
%r_ECI                  [3 x N]                         Final position
%                                                       vector in ECI
%
%v_ECI                  [3 x N]                         Final velocity
%                                                       vector in ECI
%--------------------------------------------------------------------------
% Programmed by Darin Koblick 03-04-2012                                  %
%-------------------------------------------------------------------------- 


if ~exist('mu','var'); mu = 398600.4418; end

tol = 1e-9;
v0Mag = sqrt(sum(v0.^2,1));  r0Mag = sqrt(sum(r0.^2,1));
alpha = -(v0Mag.^2)./mu + 2./r0Mag; 

%% Compute initial guess (X0) for Newton's Method
X0 = NaN(size(t));

%Check if there are any Eliptic/Circular orbits
idx = alpha > 0.000001;
if any(idx)
    X0(idx) = sqrt(mu).*t(idx).*alpha(idx); 
end

%Check if there are any Parabolic orbits
idx = abs(alpha) < 0.000001;
if any(idx)
   h = cross(r0(:,idx),v0(:,idx)); hMag = sqrt(sum(h.^2,1));
   p = (hMag.^2)./mu; s = acot(3.*sqrt(mu./(p.^3)).*t(idx))./2;
   w = atan(tan(s).^(1/3)); X0(idx) = sqrt(p).*2.*cot(2.*w);
end

%Check if there are any Hyperbolic orbits
idx = alpha < -0.000001;
if any(idx)
   a = 1./alpha(idx);
   X0(idx) = sign(t(idx)).*sqrt(-a).*...
       log(-2.*mu.*alpha(idx).*t(idx)./ ...
       (dot(r0(:,idx),v0(:,idx))+sign(t(idx)).*sqrt(-mu.*a).*...
       (1-r0Mag(idx).*alpha(idx))));
end

%% Newton's Method to converge on solution
% Declare Constants that do not need to be computed within the while loop
err = Inf;
dr0v0Smu = dot(r0,v0)./sqrt(mu);
Smut = sqrt(mu).*t;
while any(abs(err) > tol)
    X02 = X0.^2;
    X03 = X02.*X0;
    psi = X02.*alpha;
    [c2,c3] = c2c3(psi);
    X0tOmPsiC3 = X0.*(1-psi.*c3);
    X02tC2 = X02.*c2;
    r = X02tC2 + dr0v0Smu.*X0tOmPsiC3 + r0Mag.*(1-psi.*c2);
    Xn = X0 + (Smut-X03.*c3-dr0v0Smu.*X02tC2-r0Mag.*X0tOmPsiC3)./r;
    err = Xn-X0; X0 = Xn;
end
f = 1 - (Xn.^2).*c2./r0Mag; 

g = t - (Xn.^3).*c3./sqrt(mu);

gdot = 1 - c2.*(Xn.^2)./r; 

fdot = Xn.*(psi.*c3-1).*sqrt(mu)./(r.*r0Mag);

r = bsxfun(@times,f,r0) + bsxfun(@times,g,v0);
v = bsxfun(@times,fdot,r0) + bsxfun(@times,gdot,v0);

%%5 Ensure Solution Integrity
%idx = round((f.*gdot - fdot.*g)./tol).*tol ~= 1; r(:,idx) = NaN; v(:,idx) = NaN;
end

function [c2,c3] = c2c3(psi)
%Vallado pg. 71 Algorithm 1
c2 = NaN(size(psi));
c3 = NaN(size(psi));
idx = psi > 1e-6;
if any(idx)
    c2(idx) = (1-cos(sqrt(psi(idx))))./psi(idx);
    c3(idx) = (sqrt(psi(idx))-sin(sqrt(psi(idx))))./sqrt(psi(idx).^3);
end
idx = psi < -1e-6;
if any(idx)
    c2(idx) = (1 - cosh(sqrt(-psi(idx))))./psi(idx);
    c3(idx) = (sinh(sqrt(-psi(idx)))-sqrt(-psi(idx)))./sqrt(-psi(idx).^3);
end
idx = abs(psi) <= 1e-6;
if any(idx)
    c2(idx) = 0.5;
    c3(idx) = 1/6;
end
end

function [statedot] = orbitPROP(t,state,mumoon)
% All vectors must be input as COLUMNS 

% state includes the following: 
% state(1) = rAx eci
% state(2) = rAy eci
% state(3) = rAz eci
% state(4) = vAx eci
% state(5) = vAy eci
% state(6) = vAz eci
% state(7) = deltarx


R = state(1:3);
V = state(4:6);

dx = state(4);
dy = state(5);
dz = state(6);

r = norm([state(1) state(2) state(3)]);
ddx = -mumoon*state(1)/r^3;
ddy = -mumoon*state(2)/r^3;
ddz = -mumoon*state(3)/r^3;

statedot =  [dx;dy;dz;ddx;ddy;ddz];  % Propagation of R and V vectors
    

end

function [statedot] = orbitPropNonImpulsive(t,state,mumoon)%,T,Isp)
% All vectors must be input as COLUMNS 

% state includes the following: 
% state(1) = rAx eci
% state(2) = rAy eci
% state(3) = rAz eci
% state(4) = vAx eci
% state(5) = vAy eci
% state(6) = vAz eci
% state(7) = mdot

g0 = .00162;
T = 14/(1000*1000); % mN
Isp = 1100; % seconds

R = state(1:3);
V = state(4:6);

dx = state(4);
dy = state(5);
dz = state(6);

m = state(7);
v = norm([dx dy dz]);
r = norm([state(1) state(2) state(3)]);

ddx = -mumoon*state(1)/(r^3) + (T/m)*(dx/v);
ddy = -mumoon*state(2)/(r^3) + (T/m)*(dy/v);
ddz = -mumoon*state(3)/(r^3) + (T/m)*(dz/v);

mdot = -T/(Isp*g0);

statedot =  [dx;dy;dz;ddx;ddy;ddz;mdot];  % Propagation of R, V, and m

end

function [R,V] = coes2rvKAMME(h,inc,RAAN,e,per,theta,mumoon)
% h [km^2/s] Specific angular momentum
% i [rad] Inclination
% RAAN [rad] Right ascension (RA) of the ascending node
% e Eccentricity
% per [rad] Argument of perigee
% theta [rad] True anomaly
% mumoon = 4903; Moon’s gravitational parameter [km^3/s^2]

% State Vectors in Perifocal coordinates
rx = h^2/mumoon*(1/(1 + e*cos(theta)))*[cos(theta);sin(theta);0];
vx = mumoon/h*[-sin(theta); (e +cos(theta));0];

% Direction cosine matrix
DCM = [cos(per), sin(per),0;-sin(per),cos(per),0;0,0,1]*...
 [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)]*...
 [cos(RAAN), sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];

% Transformation Matrix
Dcm = inv(DCM);

% ECI R
R = Dcm*rx;

% ECI V
V = Dcm*vx;

end

function [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(r,v,mu)
%----------------------- Begin Code Sequence -----------------------------%
% Purpose:                                                                %
% Convert a given set of state vectors in ECI reference frame to orbital  %
% elements.                                                               %
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
% Outputs:
%---------                                                                %
%a                      [1 x N]                         Semi-Major Axis
%                                                       (km)
%
%eMag                   [1 x N]                         Eccentricity
%                                                       (unitless)
%
%i                      [1 x N]                         inclination
%                                                       (radians)
%
%O                      [1 x N]                         Right Ascention of
%                                                       the ascending node
%                                                       (radians)
%
%o                      [1 x N]                         Argument of perigee
%                                                       (radians)
%
%nu                      [1 x N]                         Mean Anomaly
%                                                       (radians)
%
%truLon                 [1 x N]                         True Longitude
%                                                       (radians)
%
%argLat                 [1 x N]                         Argument of Latitude
%                                                       (radians)
%
%lonPer                 [1 x N]                         Longitude of Periapse
%                                                       (radians)
%
%p                      [1 x N]                         Semilatus Rectum
%                                                       (km)
%
% References:
%-------------
%Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
%
% Function Dependencies:
%------------------
%None
%------------------------------------------------------------------       %
% Programed by Darin Koblick  03-04-2012                                  %
% Updated to address circular equatorial orbits       12/12/2013          %
%------------------------------------------------------------------       %

if ~exist('mu','var');  t = getConst(); mu = t.Earth.Mu; end
%Specific angular momentum
h = cross(r,v);
%h = h'; % mac does some whack stuff with cross
n = cross(repmat([0;0;1],[1,size(r,2)]),h);
%n = n';  % mac does some whack stuff with cross
nMag = sqrt(sum(n.^2,1));
vMag = sqrt(sum(v.^2,1)); 
rMag = sqrt(sum(r.^2,1)); 
hMag = sqrt(sum(h.^2,1));
e = (1./mu).*(bsxfun(@times,(vMag.^2 - mu./rMag),r) - bsxfun(@times,dot(r,v),v)); 
eMag = sqrt(sum(e.^2,1));
zeta = (vMag.^2)./2 - mu./rMag;

%Special Procedure when we have a parabolic orbit
idx = eMag ~= 1;
a = NaN(size(eMag));
p = NaN(size(eMag));
if any(idx)
    a(idx) = -mu./(2.*zeta(idx)); 
    p = a(idx).*(1-eMag(idx).^2); 
else
    a(idx) = Inf; 
    p(idx) = (hMag(idx).^2)./mu; 
end

%Compute the angles
i = acos(h(3,:)./hMag); 
O = acos(n(1,:)./nMag);
o = acos(dot(n,e)./(nMag.*eMag));
nu = acos(dot(e,r)./(eMag.*rMag));
lonPer = acos(e(1,:)./eMag);
argLat = acos(dot(n,r)./(nMag.*rMag));
truLon = acos(r(1,:)./rMag);

%Account for those cases where satellite is in circular orbit
         O(n(1,:) == 0) = 0;
       o(dot(n,e) == 0) = 0;
    lonPer(e(1,:) == 0) = 0;
      nu(dot(e,r) == 0) = 0;
  argLat(dot(n,r) == 0) = 0;
  
%Apply Quadrant Checks to All Determined Angles
idx = n(2,:) < 0; if any(idx);  O(idx) = 2*pi - O(idx);  end
idx = e(3,:) < 0; if any(idx); o(idx) = 2*pi - o(idx); end
idx = dot(r,v) < 0; if any(idx); nu(idx) = 2*pi - nu(idx); end
idx = e(2,:) < 0; if any(idx); lonPer(idx) = 2*pi-lonPer(idx);  end
idx = r(3,:) < 0; if any(idx); argLat(idx) = 2*pi - argLat(idx); end
idx = r(2,:) < 0; if any(idx); truLon(idx) = 2*pi - truLon(idx); end
end


function [AnonImpulsive] = nonImpulsive(state) 

g0 = .00162;
T = 14/(1000*1000); % mN
Isp = 1100; % seconds

g0 = .00162;
T = 14/(1000*1000); % mN
Isp = 1100; % seconds

R = state(1:3);
V = state(4:6);

dx = state(4);
dy = state(5);
dz = state(6);

m = state(7);
v = norm([dx dy dz]);
r = norm([state(1) state(2) state(3)]);

ddx = (T/m)*(dx/v);
ddy = (T/m)*(dy/v);
ddz = (T/m)*(dz/v);

AnonImpulsive = [ddx;ddy;ddz];

end

function coe = rv2coes(R,V,mu)
eps = 1.e-10;
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);

%...Equation 4.7:
incl = acos(H(3)/h);

%...Equation 4.8:
N = cross([0 0 1],H);
n = norm(N);
%...Equation 4.9:

if n ~= 0
RA = acos(N(1)/n);
if N(2) < 0
RA = 2*pi - RA;
end
else
RA = 0;
end

%...Equation 4.10:
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

%...Equation 4.12 (incorporating the case e = 0):
if n ~= 0
if e > eps
w = acos(dot(N,E)/n/e);
if E(3) < 0
w = 2*pi - w;
end
else
w = 0;
end
else
w = 0;
end


%...Equation 4.13a (incorporating the case e = 0):
if e > eps
    TA = acos(dot(E,R)/e/r);
if vr < 0
    TA = 2*pi - TA;
end
else
cp = cross(N,R);
if cp(3) >= 0
TA = acos(dot(N,R)/n/r);
else
TA = 2*pi - acos(dot(N,R)/n/r);
end
end

%...Equation 4.62 (a < 0 for a hyperbola):
a = h^2/mu/(1 - e^2);
coe = [h e RA incl w TA a];

end

