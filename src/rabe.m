% Implement Rabe's approach to finding trojan orbits

% Sun - Jupiter system
M = 1/1047.355;         % Mass of Jupiter (relative to Sun of mass 1)
N = sqrt(1 + M);        % Angular momentum of Jupiter
P = 2*pi/N;             % Period of Jupiter

% Trojan starting conditions
a0;
b0;
a1;
b1;

% Find other coefficients
c0 = sqrt(a0^2 + b0^2);
d0 = sqrt((a0-1)^2 + b0^2);
e0 = M * (c0^-3 - 1);
f0 = d0^-3 - 1;

a = [a0 a1];
b = [b0 b1];
c = [c0];
d = [d0];
for n=1:10
    suma = 0;
    sumb = 0;
    sumc = 0;
    sumd = 0;
    sume = 0;
    sumf = 0;
    for v=0:n
        suma = suma + a(v+1)*a(n-v+1);
        sumb = sumb + b(v+1)*b(n-v+1);
        suma = suma + (v+1)*c(n-v+1);
        suma = suma + a(v+1)*a(n-v+1);
        
    end
    c(n+1) = () / (2*c0);
end
