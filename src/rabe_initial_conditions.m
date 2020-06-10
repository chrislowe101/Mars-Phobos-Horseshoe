% Implement Section 3 of Rabe, to find approximate initial conditions (a0, b0, a1, b1) for a periodic orbit

lambda = 0.06;    % Parameter, distance from L4/L5 Lagrange point along line from that point to Sun

% Sun - Jupiter system
M = 1/1047.355;         % Mass of Jupiter (relative to Sun of mass 1)
N = sqrt(1 + M);        % Angular momentum of Jupiter
P = 2*pi/N;             % Period of Jupiter

% Initial position
d0 = 1 + lambda;
a0 = 1/2 * (1 - lambda);
b0 = sqrt(3)/2 * (1 + lambda);      % For one lagrange point only - negative is for other lagrange point (symmetry about x)

% Estimate initial velocity
c0 = sqrt(a0^2 + b0^2);
d0 = sqrt((a0-1)^2 + b0^2);
e0 = M * (c0^-3 - 1);
f0 = d0^-3 - 1;

A = (a0 - 1)*f0 + a0*e0;
B = b0*(e0 + f0);
z = -A/B;
D = (3*(1+f0))/d0^2 * ((a0-1) + b0*z)^2 + (3*(e0+M))/c0^2 * (a0 + b0*z)^2 - (e0+f0) * (1+z^2);
F = (2*N*(1+z^2)*B)/D;
G = (A^2+B^2)/D;

a1 = min(abs(roots([1, F, G])));    % Two roots obtained - smaller one applicable in this case
b1 = z*a1;

% TEST - ensure output matches Table I
d0
a0
b0
a1
b1
