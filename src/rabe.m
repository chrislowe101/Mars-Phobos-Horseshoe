% Implement Rabe's approach to finding trojan orbits

% Sun - Jupiter system
M = 1/1047.355;         % Mass of Jupiter (relative to Sun of mass 1)
N = sqrt(1 + M);        % Angular momentum of Jupiter
P = 2*pi/N;             % Period of Jupiter

% Trojan starting conditions
a0 = 0.47;
b0 = 0.917986928;
a1 = 0.077336855;
b1 = 0.044371494;

% Find other initial coefficients
c0 = sqrt(a0^2 + b0^2);
d0 = sqrt((a0-1)^2 + b0^2);
e0 = M * (c0^-3 - 1);
f0 = d0^-3 - 1;

a = [a0 a1];
b = [b0 b1];
c = [c0];
d = [d0];

% Recursively find higher order coefficients
for n = 1:10            % hard limit for now, but should a stopping condtion for when all six coefficients are < some epsilon

    % Find next c coefficient
    suma = 0;
    sumb = 0;
    sumc = 0;
    for v = 0:n
        suma = suma + a(v+1)*a(n-v+1);
        sumb = sumb + b(v+1)*b(n-v+1);
    end
    if (n-1) >= 1
        for v = 1:(n-1)
            sumc = sumc + c(v+1)*c(n-v+1);
        end
    end
    c(n+1) = (suma + sumb - sumc) / (2*c(1));
    
    % Find next d coefficient
    sumc = 0;
    sumd = 0;
    for v = 0:n
        sumc = sumc + c(v+1)*c(n-v+1);
    end
    for v = 0:(n-1)
        sumd = sumd + d(v+1)*d(n-v+1);
    end
    d(n+1) = (sumc - sumd - 2*a(n+1)) / (2*d(1));
    
    % Find next e coefficient
    sume1 = 0;
    sume2 = 0;
    for v = 1:n
        sume1 = sume1 + v*c(v+1)*e(n-v+1);
    end
    if (n-1) >= 1
        for v = 1:(n-1)
            sume2 = sume2 + v*e(v+1)*c(n-v+1);
        end
    end
    e(n+1) = (3*sume1 + sume2 + 3*M*n*c(n+1)) / (-n*c(1));
    
    % Find next f coefficient
    sumf1 = 0;
    sumf2 = 0;
    for v = 1:n
        sumf1 = sumf1 + v*d(v+1)*f(n-v+1);
    end
    if (n-1) >= 1
        for v = 1:(n-1)
            sumf2 = sumf2 + v*f(v+1)*d(n-v+1);
        end
    end
    f(n+1) = (3*sumf1 + sumf2 + 3*n*d(n+1)) / (-n*d(1));
    
    % Find next a coefficient
    suma = 0;
    for v = 0:(n-1)
        suma = suma + a(v+1)*(e(n-v) + f(n-v));
    end
    a(n+2) = (suma - 2*N*n*b(n+1) - f(n)) / (-n*(n+1));
    
    % Find next b coefficient
    sumb = 0;
    for v = 0:(n-1)
        sumb = sumb + b(v+1)*(e(n-v) + f(n-v));
    end
    b(n+2) = (sumb + 2*N*n*a(n+1)) / (-n*(n+1));
    
end

% TEST - compare to paper results (starting conditions should yield d0 = 1.06 and residuals as per Table II, Table VIII
a
b

% Find position and velocity at this time step
p = sum(a);
q = sum(b);
p_dot = 0;
q_dot = 0;
for n = 1:length(a)
    p_dot = p_dot + n*a(n+1);
    q_dot = q_dot + n*b(n+1);
end

% Calculate Jacobi constant
r = sqrt(p^2 + q^2);        % Distance from Jupiter to Trojan
s = sqrt((p-1)^2 + q^2);    % Distance from Sun to Trojan
C = M*(r^2 + 2/r) + (s^2 + 2/s) - p_dot^2 - q_dot^2;
