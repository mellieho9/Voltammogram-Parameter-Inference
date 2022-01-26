%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CVsim.m - Cyclic voltammetry simulation
% Meliora Ho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INDEPENDENT VARIABLES %%
C      = 1.0;    % [=] mol/L, initial concentration of O. Default = 1.0
D      = 1E-5;   % [=] cm^2/s, O & R diffusion coefficient. Default = 1E-5
etai   = +0.2;   % [=] V, initial overpotential (relative to redox potential). Default = +0.2
etaf   = -0.2;   % [=] V, final overpotential (relative to redox potential). Default = -0.2
v      = 1E-3;   % [=] V/s, sweep rate. Default = 1E-3
n      = 1.0;    % [=] number of electrons transfered. Default = 1
alpha  = 0.5;    % [=] dimensionless charge-transfer coefficient. Default = 0.5
k0     = 1E-2;   % [=] cm/s, electrochemical rate constant. Default = 1E-2
kc     = 1E-3;   % [=] 1/s, chemical rate constant. Default = 1E-3
T      = 298.15; % [=] K, temperature. Default = 298.15

%% PHYSICAL CONSTANTS %%
F      = 96485;   % [=] C/mol, Faraday's constant
R      = 8.3145;  % [=] J/mol-K, ideal gas constant
f      = F/(R*T); % [=] 1/V, normalized Faraday's constant at room temperature

%% SIMULATION VARIABLES %%
L      = 500;    % [=] number of iterations per t_k. Default = 500
DM     = 0.45;   % [=] model diffusion coefficient. Default = 0.45

%% DERIVED CONSTANTS %%
tk  = 2*(etai-etaf)/v;    % [=] s, characteristic exp. time. In this case, total time of fwd and rev scans
Dt  = tk/L;               % [=] s, delta time
Dx  = sqrt(D*Dt/DM);      % [=] cm, delta x 
j   = ceil(4.2*L^0.5)+5;  % number of boxes

%% PRE-INITIALIZATION %%
C = C / 1000;           % Convert C from mol/L to mol/cm3
k = 0:L;                % time index vector
t = Dt * k;             % time vector
eta1 = etai - v*t;      % overpotential vector, negative scan
eta2 = etaf + v*t;      % overpotential vector, positive scan
eta = [eta1(eta1>etaf) eta2(eta2<=etai)]'; % overpotential scan, both directions
Enorm = eta*f;          % normalized overpotential
kf = k0.*exp(  -alpha *n*Enorm); % [=] cm/s, fwd rate constant
kb = k0.*exp((1-alpha)*n*Enorm); % [=] cm/s, rev rate constant

O = C*ones(L+1,j); % [=] mol/cm^3, concentration of O
R = zeros(L+1,j);  % [=] mol/cm^3, concentration of R
JO = zeros(1,L+1); % [=] mol/cm^2-s, flux of O at the surface

%% START SIMULATION %%
% i1 = time index. i2 = distance index
for i1 = 1:L
    % Update bulk concentrations of O and R
    for i2 = 2:j-1
        O(i1+1,i2) = O(i1,i2) + DM*(O(i1,i2+1)+O(i1,i2-1)-2*O(i1,i2));

        R(i1+1,i2) = R(i1,i2) + DM*(R(i1,i2+1)+R(i1,i2-1)-2*R(i1,i2)) ...
            - km * R(i1,i2);
    end

    % Update flux
    JO(i1+1)   = ( kf(i1+1).*O(i1+1,2) - kb(i1+1).*R(i1+1,2) ) ./ (1 + Dx/D*(kf(i1+1) + kb(i1+1)) );

    % Update surface concentrations
    O(i1+1,1) = O(i1+1,2) - JO(i1+1)*(Dx/D);
    R(i1+1,1) = R(i1+1,2) + JO(i1+1)*(Dx/D) - km*R(i1+1,1);
end

% Calculate current density, Z, from flux of O
Z = -n.*F.*JO * 1000; % [=] A/cm^2 -> mA/cm^2, current density

%% PLOT RESULTS %%
if length(eta) > length(Z)
    eta = eta(1:end-1);
end

plot(eta,Z)
xlabel('Overpotential (V)'), ylabel('Current density (mA/cm^2)')
