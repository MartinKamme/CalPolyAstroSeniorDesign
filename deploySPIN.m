%% AERO 449 Martin Kamme Deployment With Spin
close all
clear all
clc

mumoon = 4903;                              % Moon's Gravitational Parameter km^3/s^2
T = 14/(1000*1000);                         % Thrust [kg*km/s^2]
Isp = 1100;                                 % Isp [seconds]
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


% LCV Orbit----------------------------------------------------------------
below_op = 120;                             % LCV alt below Operational orbit [km]

inc_lcv = deg2rad(55);                      % Inclination [inc]
ecc_lcv = 0;                                % Eccentricity
RAAN_lcv = 0;                               % Right Ascension of Ascending Node [rad]
r_lcv = a_op - below_op;                    % Orbital Radius [km]
T_lcv = ((2*pi)/sqrt(mumoon))*(r_lcv^(3/2));% Period [seconds]
V_lcv = sqrt(mumoon/r_lcv);                 % Velocity [km/s]
h_lcv = sqrt(r_lcv*mumoon*(1-ecc_lcv^2));   % Angular Momentum [km^2/s]
theta_lcv = deg2rad(0);                     % True anomaly [rad]
per_lcv = deg2rad(0);                       % Argument of perilune [rad] define as being located 0 degrees from the ascending node
n_lcv = 2*pi/(T_lcv);                       % Mean Motion [rad/s] 
thetaINT = deg2rad(270);

LCVplotRadius = r_lcv;

[R_lcv,V_lcv] = coes2rvKAMME(h_lcv,inc_lcv,RAAN_lcv,ecc_lcv,per_lcv,theta_lcv,mumoon);
[R_lcvRSW,V_lcvRSW] = ECI2RSW(R_lcv,V_lcv);


% State Vectors and Things-------------------------------------------------
[R_lss,V_lss] = coes2rvKAMME(h_op,inc_op,RAAN_op,ecc_op,per_op,thetaINT,mumoon);


% TAKE INTO ACCOUNT THE SPIN AND SPRINGINESS ------------------------------

% Spin
rpm = 1; % rotation per minute
LCVradius = 31/(100*2.54);                  % cm
angVel = 2*pi*rpm/60;                       % rad/s
releaseVel = (LCVradius*angVel)/1000;       % km/s

% Springiness
e = 0.85 ;                                  % J
numSprings = 6 ;
massLSS = 179 ;                             % kg
massLCV = 3896 ;                            % kg
vrel = sqrt((numSprings*(massLSS+massLCV)*2*0.85)/(massLCV*massLSS)); % m/s
vLSS_separation = (vrel/(1-(massLSS/massLCV)))/1000;                  % km/s

% Combination of spin and release
vUnit = [releaseVel ;0 ;vLSS_separation];  % R S W
angles = linspace(0,2*pi);
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(R_lcv,V_lcv,mumoon);


for ii = 1:length(angles)
    velRELEASE(:,ii) = rotateVrelease(vUnit,angles(ii));


V_lcvSPIN = V_lcvRSW + velRELEASE(:,ii);

[R_SPIN_ECI,V_SPIN_ECI] = RSW2ECI(R_lcvRSW,V_lcvSPIN,RAAN_lcv,inc_lcv,argLat);

% Propagate non impulsive with Encke's
dr = [0;0;0];
dv = [0;0;0];
dt = 10; % set time step as 10 secondS

rpert = R_SPIN_ECI;
vpert = V_SPIN_ECI;
rosc = R_SPIN_ECI;
vosc = V_SPIN_ECI;

radius = 12;
t = 0; % initiliaze time
i = 1; % counter for tracking values
m = 179; 

% Propagate the LSS WITH SPIN
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
    rplotLSSSPIN(:,i) = rpert.*1;
    vplotLSSSPIN(:,i) = vpert.*1;
    coe = rv2coes(rplotLSSSPIN(:,i),vplotLSSSPIN(:,i),mumoon);
    radius = coe(1,7);
    
    % rectify every step
    rosc = rpert;
    vosc = vpert; 
    dr = [0;0;0];
    dv = [0;0;0];
    
    i = i+1;
    
end



% figure
% hold on 
% plot3(rplotLSS(1,:),rplotLSS(2,:),rplotLSS(3,:))
% plot3(rplotLSSSPIN(1,:),rplotLSSSPIN(2,:),rplotLSSSPIN(3,:))
% legend('LSS','LSS with spin')


% Important Post Calculations
coesSPIN(ii,:) = coe;

end

% Propagate non impulsive with Encke's
dr = [0;0;0];
dv = [0;0;0];
dt = 10; % set time step as 10 seconds

rpert = R_lcv;
vpert = V_lcv;
rosc = R_lcv;
vosc = V_lcv;

radius = 12;
t = 0; % initiliaze time
i = 1; % counter for tracking values
m = 179; 

% Propagate the LSS WITHOUT SPIN
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
    rplotLSS(:,i) = rpert.*1;
    vplotLSS(:,i) = vpert.*1;
    coe = rv2coes(rplotLSS(:,i),vplotLSS(:,i),mumoon);
    radius = coe(1,7);

    % rectify every step
    rosc = rpert;
    vosc = vpert; 
    dr = [0;0;0];
    dv = [0;0;0];
    
    i = i+1;
    
end
coeNOSPIN = coe;

%%
close all 
clc

figure
title('Impact of LCV Spin on Angular Momentum')
hold on 
grid on
plot(200,coeNOSPIN(1),'r*')
plot(rad2deg(angles),coesSPIN(:,1))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('Momentum')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on Eccentricity')
hold on 
grid on
plot(200,coeNOSPIN(2),'r*')
plot(rad2deg(angles),coesSPIN(:,2))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('Ecentrictity')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on RAAN')
hold on 
grid on
plot(200,rad2deg(coeNOSPIN(3)),'r*')
plot(rad2deg(angles),rad2deg(coesSPIN(:,3)),'bo')
xlabel('Release Velocity Angle from R (degrees)')
ylabel('RAAN (degrees)')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on Inclination')
hold on 
grid on
plot(200,rad2deg(coeNOSPIN(4)),'r*')
plot(rad2deg(angles),rad2deg(coesSPIN(:,4)))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('Inclination (degrees)')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on Argument of Periselene')
hold on 
grid on
plot(200,rad2deg(coeNOSPIN(5)),'r*')
plot(rad2deg(angles),rad2deg(coesSPIN(:,5)))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('Argument of Periselene (degrees)')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on True Anomaly')
hold on 
grid on
plot(200,rad2deg(coeNOSPIN(6)),'r*')
plot(rad2deg(angles),rad2deg(coesSPIN(:,6)))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('True Anomaly (degrees)')
legend('No spin, (not correlated to 200 degrees)')

figure
title('Impact of LCV Spin on Semi Major Axis')
hold on 
grid on
plot(200,coeNOSPIN(7),'r*')
plot(rad2deg(angles),coesSPIN(:,7))
xlabel('Release Velocity Angle from R (degrees)')
ylabel('Semi Major Axis')
legend('No spin, (not correlated to 200 degrees)')

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


% Vallado Hill Stuff
function [rHill,vHill] = ECI2Hill_Vectorized(rTgt,vTgt,rChase,vChase)
% Purpose:
% Convert those position (ECI) and velocity (ECI) into Hill's reference
% frame using both the target and the chaser position/velocity data
%
% Inputs:
%rTgt                       [3 x N]                 ECI Position vector of
%                                                   reference frame (km)
%
%vTgt                       [3 x N]                 ECI Velocity vector of
%                                                   reference frame (km/s)
%rChase                     [3 x N]
%
%vChase                     [3 x N]
%
% Outputs:
%rHill                      [3 x N]                 Hill's relative
%                                                   position vector (km)
%
%vHill                      [3 x N]                 Hill's relative
%                                                   velocity vector (km/s)
% References:
% Vallado 2007.
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
%Declare Local Functions
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
rTgtMag = sqrt(sum(rTgt.^2,1));
rChaseMag = sqrt(sum(rChase.^2,1));
vTgtMag = sqrt(sum(vTgt.^2,1));
%Determine the RSW transformation matrix from the target frame of reference
RSW = ECI2RSW(rTgt,vTgt);
%Use RSW rotation matrix to convert rChase and vChase to RSW
r_Chase_RSW = matrixMultiply(RSW,rChase);
v_Chase_RSW = matrixMultiply(RSW,vChase);
%Find Rotation angles to go from target to interceptor
phi_chase = asin(r_Chase_RSW(3,:)./rChaseMag);
lambda_chase = atan2(r_Chase_RSW(2,:),r_Chase_RSW(1,:));
CPC = cos(phi_chase);     SPC = sin(phi_chase);
SLC = sin(lambda_chase);  CLC = cos(lambda_chase);
%Find Position component rotations
rHill = cat(1,rChaseMag-rTgtMag, ...
              lambda_chase.*rTgtMag, ...
              phi_chase.*rTgtMag);
%Find the rotation matrix RSW->SEZ of chaser
RSW_SEZ = zeros(3,3,size(rTgtMag,2));
RSW_SEZ(1,1,:) = SPC.*CLC;  RSW_SEZ(1,2,:) = SPC.*SLC;  RSW_SEZ(1,3,:) = -CPC;
RSW_SEZ(2,1,:) = -SLC;  RSW_SEZ(2,2,:) = CLC;
RSW_SEZ(3,1,:) = CPC.*CLC;  RSW_SEZ(3,2,:) = CPC.*SLC;  RSW_SEZ(3,3,:) = SPC;
%Find the velocity component of positions using the angular rates in SEZ frame
v_Chase_SEZ = matrixMultiply(RSW_SEZ,v_Chase_RSW);
vHill = cat(1,v_Chase_SEZ(3,:), ...
              rTgtMag.*(v_Chase_SEZ(2,:)./(rChaseMag.*CPC)-vTgtMag./rTgtMag), ...
              -rTgtMag.*v_Chase_SEZ(1,:)./rChaseMag);
end

function [rRSW,vRSW] = ECI2RSW(rECI,vECI)
% Purpose:
%Convert ECI Coordinates to RSW Coordinates, also, return the
%transformation matrix T in which to take a given set of coordinates in ECI
%and convert them using the same RSW reference frame.
%
% Inputs:
% rECI              [3 x N]                     ECI position Coordinates in
%                                               km
%
% vECI              [3 x N]                     ECI velocity Coordinates in
%                                               km/s
%
% Outputs:
% T                 [3 x 3 x N]                 Transformation matrix
%                                               necessary to go from
%                                               rECI -> rRSW
%
% rRSW              [3 x N]                     RSW Position Coordinates
%                                               km
%
% vRSW              [3 x N]                     RSW Velocity Coordinates
%                                               km/s
%
% References:
% Vallado pg. 173
% Vallado rv2rsw.m code dated 06/09/2002
%
%Programmed by: Darin C Koblick                 11/29/2012
% Begin Code Sequence
if nargin == 0
    rECI =  repmat([6968.1363,1,2]',[1 10]);
    vECI =  repmat([3,7.90536615282099,4]',[1 10]);
    [T,rRSW,vRSW] = ECI2RSW(rECI,vECI);
    return;
end
%Declared Internal function set:
unitv = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2,1)));
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
%Find the Radial component of the RIC position vector
rvec = unitv(rECI);
%Find the cross-track component of the RIC position vector
wvec = unitv(cross(rECI,vECI));
%Find the along-track component of the RIC position vector
svec = unitv(cross(wvec,rvec));
%Create the transformation matrix from ECI to RSW
T = NaN(3,3,size(rECI,2));
T(1,1,:) = rvec(1,:); T(1,2,:) = rvec(2,:); T(1,3,:) = rvec(3,:);
T(2,1,:) = svec(1,:); T(2,2,:) = svec(2,:); T(2,3,:) = svec(3,:);
T(3,1,:) = wvec(1,:); T(3,2,:) = wvec(2,:); T(3,3,:) = wvec(3,:);
%Find the position and velocity vectors in the RSW reference frame!
rRSW = matrixMultiply(T,rECI);
vRSW = matrixMultiply(T,vECI);
end

function [rHill,vHill] = CWHPropagator(rHillInit,vHillInit,omega,t)
% Purpose:
% Take initial position and velocity coordinates in the Hill reference frame
% and propagate them using the Clohessy-Wiltshire Hill Linearize equation
% of motion.
%
% Inputs:
%rHillInit                  [3 x 1]                 Hill Position vector
%                                                   (km) / (m)
%
%vHillInit                  [3 x 1]                 Hill Velocity vector of
%                                                   (km/s) / (m/s)
%
%omega                       double                 Orbital Angular Rate
%                                                   of the target
%                                                   (rad/s)
%                                                   Should be close to
%                                                   circular for linear propagation
%                                                   error to be low.
%
%t                          [1 x N]                 Propagation Time in
%                                                   seconds
%                                                   
%
%
%
% Outputs:
%rHill                       [3 x N]                Propagated Hill
%                                                   Position vector (km) /
%                                                   (m/s)
%
%vHill                       [3 x N]                Propagated Hill
%                                                   Velocity vector (km/s)
%                                                   / (m/s)
%
%
% References:
%
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
x0 = rHillInit(1,:); y0 = rHillInit(2,:); z0 = rHillInit(3,:);
x0dot = vHillInit(1,:); y0dot = vHillInit(2,:); z0dot = vHillInit(3,:);
rHill = [(x0dot./omega).*sin(omega.*t)-(3.*x0+2.*y0dot./omega).*cos(omega.*t)+(4.*x0+2.*y0dot./omega)
        (6.*x0+4.*y0dot./omega).*sin(omega.*t)+2.*(x0dot./omega).*cos(omega.*t)-(6.*omega.*x0+3.*y0dot).*t+(y0-2.*x0dot./omega)
        z0.*cos(omega.*t)+(z0dot./omega).*sin(omega.*t)];
vHill = [x0dot.*cos(omega.*t)+(3.*omega.*x0+2.*y0dot).*sin(omega.*t)
        (6.*omega.*x0 + 4.*y0dot).*cos(omega.*t) - 2.*x0dot.*sin(omega.*t)-(6.*omega.*x0 + 3.*y0dot)
        -z0.*omega.*sin(omega.*t)+z0dot.*cos(omega.*t)];
end

function [rInt,vInt] = Hill2ECI_Vectorized(rTgt,vTgt,rHill,vHill)
% Purpose:
% Convert those position (rHill) and velocity (vHill) values back into an
% ECI coordinate frame of reference using the reference satellite
% (rTgt,vTgt) position and velocity data.
%
% Inputs:
%rTgt                       [3 x N]                 ECI Position vector of
%                                                   reference frame (km)
%
%vTgt                       [3 x N]                 ECI Velocity vector of
%                                                   reference frame (km/s)
%
%rHill                      [3 x N]                 Hill's relative
%                                                   position vector (km)
%
%vHill                      [3 x N]                 Hill's relative
%                                                   velocity vector (km/s)
%
%
%
% Outputs:
%rInt                       [3 x N]
%
%vInt                       [3 x N]
%
%
% References:
% Vallado 2007.
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
%Declare Local Functions
rTgtMag = sqrt(sum(rTgt.^2,1));
vTgtMag = sqrt(sum(vTgt.^2,1));
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
%Find the RSW matrix from the target ECI positions
RSW = ECI2RSW(rTgt,vTgt); rIntMag = rTgtMag + rHill(1,:);
%Compute rotation angles to go from tgt to interceptor
lambda_int = rHill(2,:)./rTgtMag; phi_int = sin(rHill(3,:)./rTgtMag);
CLI = cos(lambda_int); SLI = sin(lambda_int); CPI = cos(phi_int); SPI = sin(phi_int);
%find rotation matrix to go from rsw to SEZ of inerceptor
RSW_SEZ = zeros(3,3,size(rTgt,2));
RSW_SEZ(1,1,:) = SPI.*CLI;  RSW_SEZ(1,2,:) = SPI.*SLI; RSW_SEZ(1,3,:) = -CPI;
RSW_SEZ(2,1,:) = -SLI;      RSW_SEZ(2,2,:) = CLI;      RSW_SEZ(3,1,:) = CPI.*CLI;
                            RSW_SEZ(3,2,:) = CPI.*SLI; RSW_SEZ(3,3,:) = SPI;
%Find velocity component positions by using angular rates in SEZ frame
vIntSEZ = cat(1,-rIntMag.*vHill(3,:)./rTgtMag, ...
                 rIntMag.*(vHill(2,:)./rTgtMag + vTgtMag./rTgtMag).*CPI, ...
                 vHill(1,:));
vInt = matrixMultiply(permute(RSW,[2 1 3]), ...
       matrixMultiply(permute(RSW_SEZ,[2 1 3]), ...
       vIntSEZ));
%Find the position components
rIntRSW = bsxfun(@times,rIntMag,cat(1,CPI.*CLI, ...
                                      CPI.*SLI, ...
                                      SPI));
rInt = matrixMultiply(permute(RSW,[2 1 3]),rIntRSW);      
end


% Customized for Spin Release
function [vRELEASE] = rotateVrelease(vUnit,angle)
% Rotates a vector vUnit through an angle about the tangential velocity
% vector in RSW
% Assumes that S is the y axis 

% Input:

% vUnit            3x1 [km/s]
% angle            radians


T = [cos(angle) 0 -sin(angle);
     0 1 0;
     sin(angle) 0 cos(angle)];
 
 vRELEASE = T*vUnit;
end

function [rECI,vECI] = RSW2ECI(rRSW,vRSW,RAAN,inc,arglat)
% Converts r and v from RSW to ECI
% from p.169 in VALLADO

% Inputs:

rotinc = ROT1(-inc);
rotRAAN = ROT3(-RAAN);
rotLAT = ROT3(-arglat);

rECI = rotRAAN*rotinc*rotLAT*rRSW;
vECI = rotRAAN*rotinc*rotLAT*vRSW;


end

function [rot1] = ROT1(angle)
% Generates ROTATION 1 matrix

rot1 = [1 0 0;
        0 cos(angle) sin(angle)
        0 -sin(angle) cos(angle)];
end

function [rot3] = ROT3(angle)
% Generates ROTATION 3 matrix

rot3 = [cos(angle) sin(angle) 0;
        -sin(angle) cos(angle) 0;
        0 0 1];
end


