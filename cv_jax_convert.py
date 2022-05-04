## Import necessary libraries
import math
import jax.numpy as np
import matplotlib.pyplot as plt

## Independent variables
C = 1.0
D = 1e-5
etai = 0.2
etaf = 0.2
v = 1e-3
n = 1.0
alpha = 0.5
k0 = 1e-2
kc = 1e-3
T = 298.15

## Physical constants
F = 96485
R = 8.3145
f = F/(R*T)

## Simulation variables
L = 500
DM = 0.45

## Derived constants
tk = 2*(etai-etaf)/v
Dt = tk/L
Dx = math.sqrt(D*Dt/DM)
j = math.ceil(4.2*math.pow(L,0.5)) + 5

## Pre-initialization
C /= 1000
lst = array()
for x in range(L):
    lst.append(x)
k = np.array(lst)
t = Dt * k
eta1 = etai - v*t
eta2 = etaf + v*t
# Placeholder for eta = [eta1(eta1>etaf) eta2(eta2<=etai)]
Enorm = eta*f
kf = k0*math.exp(-alpha*n*Enorm)
kb = k0*math.exp((1-alpha)*n*Enorm)
0 = array()
for y in range(L+1,j):
    0[y]=C
R = array()
for z in range (L+1,j):
    R[z]=0
J0 = array()
for i in range (1,L+1):
    J0[i]=0

## Start simulation
##
	for i1 in range(1,L):
        % Update bulk concentrations of O and R
        for i2 in range(2,j-1):
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
## 

# Calculate current density, Z, from flux of O
Z = -n*F*JO*1000

# Plot results
if len(eta) > len(Z):
    eta = eta(1:end-1);
end

plt.plot(eta,Z)
plt.xlabel('Overpotential (V)')
plt.ylabel('Current density (mA/cm^2)')

