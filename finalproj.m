close all
clear all
clc

% Difine impact speed and material properties
vimp = 300;
rho = 2.7e3;
Cl = 6.42e3;
a1 = 0.06;
a2 = 2;
YM = 68e9;
nv = 0.32;

% Difine the domain -- thickness of impactor(li) and target(L-li)
l0 = 0;
li = 1.1e-3;        % Impactor/target interface
L = 2e-3;

% Difine simulation parameters
T = 6e-7;
h = 1e-6;
k = 1e-11;
vnodes = (l0:2*h:L)';          % Elements represent locations where v is stored
stnodes = (l0+h:2*h:L-h)';     % Elements represent locations where all other variables
                               % are stored

% Initial conditions for all state variables and v
P = zeros(length(stnodes),1);
F = ones(length(stnodes),1);
E = (F.*F - 1)/2;
S = zeros(length(stnodes),1);
v = zeros(length(vnodes),1);
Volexp = zeros(length(stnodes),1);

for M = 1:length(vnodes)
    if vnodes(M) <= li
        v(M) = vimp;
    else
        v(M) = 0;
    end
end

n = 1;
while (n-1) * k < T
    
    for m = 1:M
        
        % Left boundary of domain, where all (m-1/2) terms are zero
        if m == 1
            if Volexp(m,n) < 0
            qlb = rho * (Cl * a1 * abs(v(m+1,n) - v(m,n)) + a2 * (v(m+1,n) - v(m,n))^2);
            v(m,n+1) = v(m,n) + k/(h*rho) * (P(m,n) - qlb);
            else
                v(m,n+1) = v(m,n) + k/(h*rho) * (P(m,n));
            end
            
            % Right boundary of domain, where all (m+1/2) terms are zero
        else if m == M;
                if Volexp(m-1,n) < 0
                qrb = rho * (Cl * a1 * abs(v(m,n) - v(m-1,n)) + a2 * (v(m,n) - v(m-1,n))^2);
                v(m,n+1) = v(m,n) + k/(h*rho) * (-P(m-1,n) + qrb);
                else
                 v(m,n+1) = v(m,n) + k/(h*rho) * (-P(m-1,n));
                end
                % Interior
            else
                if Volexp(m,n) + Volexp(m-1,n) < 0
                ql = rho * (Cl * a1 * abs(v(m,n) - v(m-1,n)) + a2 * (v(m,n) - v(m-1,n))^2);
                qr = rho * (Cl * a1 * abs(v(m+1,n) - v(m,n)) + a2 * (v(m+1,n) - v(m,n))^2);
                v(m,n+1) = v(m,n) + k/(2*h*rho) * (P(m,n)-P(m-1,n) - qr + ql);
                else
                   v(m,n+1) = v(m,n) + k/(2*h*rho) * (P(m,n)-P(m-1,n));
                end
            end
        end
    end
    
    for l = 1:M-1
        F(l,n+1) = F(l,n) + (k/(2*h)) * (v(l+1,n+1) - v(l,n+1));
    end
    
        E(:,n+1) = (F(:,n+1) .* F(:,n+1)-1) / 2;
        S(:,n+1) = (YM / (1+nv) + YM * nv / ((1+nv)*(1-2*nv))) * E(:,n+1);
        P(:,n+1) = S(:,n+1) .* F(:,n+1);
        Volexp(:,n+1) = E(:,n+1) - E(:,n);
        n = n + 1;
end

[r, c] = size(v);

hold on
% Particles inside the material
plot((1:c)*1e-4,v(floor(r*0.75),1:c),'r','linewidth',1.5);
% Particles at free surface
plot((1:c)*1e-4,v(r,1:c),'k','linewidth',1.5);
xlabel('Time (us)','fontsize',12)
ylabel('Velocity (m/s)','fontsize',12)
ylim([-400,500])
title('Particle velocity vs. time','fontsize',12)







