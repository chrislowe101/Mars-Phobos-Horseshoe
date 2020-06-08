%% Function: [p,a,e,i,Omega,omega,nu] = RV2COE(r,v,mu)
%
% Performs conversion from position and velocity vectors to classical
% orbital elements. Includes handling of special cases (any orbit that is
% equatorial and/or circular). In these cases, the value of nu (true
% anomoly) is undefined in the normal sense and is instead assigned the
% following values:
%   For elliptical equatorial - true longitude of periapsis
%   For circular equatorial - argument of latitude
%   For circular inclined - true longitude
% 
% Inputs:       Description                             Range / Units
%   r           - Position vector (IJK frame)           [km, km, km]
%   v           - Velocity vector (IJK frame)           [km/s, km/s, km/s]
%   mu          - Gravitational parameter               km^3 / s^2
%
% Outputs:
%   p           - Semiparameter                         km
%   a           - Semimajor Axis                        km
%   e           - Eccentricity                          
%   i           - Inclination                           rad
%   Omega       - Right Ascension of the Ascending Node rad
%   omega       - Argument of Periapsis                 rad
%   nu          - True anomoly                          rad
%
% Dependencies:
%   Calls trigonometric functions
% 
% Author:       Christopher Lowe
% Revision:     21 October 2019
%
% Reference:
%   Algorithm 9
%   Vallado, Fundamentals of Astrodynamics, Fourth Edition, Page 113.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,a,e,i,Omega,omega,nu] = RV2COE(r,v,mu)

epsilon = 10^-10;

r_mag = norm(r);
v_mag = norm(v);

% Calculate vectors
h = cross(r,v);         % Angular momentum
h_mag = norm(h);
n = cross([0,0,1],h);   % Line of nodes
n_mag =norm(n);
e_vec = ((v_mag^2 - mu/r_mag)*r - dot(r,v)*v)/ mu;   % Eccentricity vector (to periapsis)
e_mag = norm(e_vec);

% Calculate specific energy
xi = 0.5*v_mag^2 - mu/r_mag;

% Calculate orbital elements for most common cases
e = e_mag;                  % Eccentricity
i = acos(h(3)/h_mag);       % Inclination

if abs(e-1) > epsilon
    a = -mu/(2*xi);         % Semimajor axis
    p = a*(1-e_mag^2);      % Semiparameter
else
    a = inf;                % Semimajor axis (parabolic case)
    p = h_mag^2/mu;         % Semiparameter (parabolic case)
end

Omega = acos(n(1)/n_mag);   % RAAN
if n(2) < 0
    Omega = 2*pi - Omega;
end

omega = acos(dot(n,e_vec)/(n_mag*e_mag));   % Argument of periapsis
if e_vec(3) < 0
    omega = 2*pi - omega;
end

nu = acos(dot(e_vec,r)/(e_mag*r_mag));      % True anomoly
if dot(r,v) < 0
    nu = 2*pi - nu;
end

% Calculate orbital elements for special cases

% Elliptical equatorial orbit (e>0,i=0) - nu represents omega_true
if (e_mag>epsilon) && (i<epsilon)
    nu = acos(e_vec(1)/e_mag);              % True longitude of periapsis
    if e_vec(2) < 0
        nu = 2*pi - nu;
    end

% Circular inclined orbit (e=0,i>0) - nu represents u
elseif (e_mag<epsilon) && (i>epsilon)
    nu = acos(dot(n,r)/(n_mag*r_mag));      % Argument of latitude
    if r(3) < 0
        nu = 2*pi - nu;
    end

% Circular equatorial orbit (e=0,i=0) - nu represents lambda_true
elseif (e_mag<epsilon) && (i<epsilon)
    nu = acos(r(1)/r_mag);                  % True longitude
    if r(2) < 0
        nu = 2*pi - nu;
    end
end
