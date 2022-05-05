## Import necessary libraries
import math
import jax.numpy as np
import matplotlib.pyplot as plt

def cv_simulator(C,D,ni,nf,v,alpha,k0,kc,T):
    ## CONSTANTS
    n = 1
    F = 96485
    R = 8.3145
    f = F/(R*T)
    ## SIMULATION VARIABLES
    L = 500
    DM = 0.45
    ## DERIVED CONSTANTS
    tk = 2*(etai-etaf)/v
    Dt = tk/L
    Dx = math.sqrt(D*Dt/DM)
    j = math.ceil(4.2*math.sqrt(L)) + 5
    ## REVERSIBILITY PARAMETERS
    ktk = kc*tk
    km = ktk/L
    Lambda = k0/math.sqrt(D*f*v)
    ## WARNINGS
    warning = ""
    if km > 0.1:
      warning += ('k_c*t_k/L equals '+ km.toString() +', ' + 'which exceeds the upper limit of 0.1 (see B&F, pg 797). Try lowering k_c\n')
    if C < 0:
    warning += 'Concentration cannot be negative\n'
    if D < 0:
      warning += 'Diffusion coefficient cannot be negative\n'
    if etai < etaf
      warning += 'Initial-final overpotential cannot be negative\n'
    if alpha < 0 || alpha > 1:
      warning += 'Alpha must range between 0 and 1\n'
    if k0 < 0:
      warning += 'Electrochemical rate constant cannot be negative\n'
    if kc < 0:
      warning += 'Chemical rate constant cannot be negative\n'
    if warning === "": 
      warning = 'No warnings'
    print(warning)
    ## PRE-INITIALIZATIONS
    k = math.range(0,L+1)
    t = math.eval('Dt .* k',{Dt: Dt, k: k})
    eta1 = math.eval('etai - v.*t',{etai: etai, v: v, t:t})
    eta2 = math.eval('etaf + v.*t',{etaf: etaf, v: v, t:t})
    ## OVERPOTENTIAL SCAN, BOTH DIRECTIONS
    eta = []
    i = 0
    curr_eta = etai
    while curr_eta > etaf:
      eta[i] = curr_eta
      i += 1
      curr_eta = eta1.loc(math_index(k))
    eta[i] = curr_eta

    C = C / 1000 ## Convert from mol/K to mol/cm^3
    Enorm = math.eval('eta.*f',{eta: eta, f: f})
    kf = math.eval('k0*exp(  -alpha *n*Enorm)',{k0:k0, alpha:alpha, n:n, Enorm:Enorm})
    kb = math.eval('k0*exp((1-alpha)*n*Enorm)',{k0:k0, alpha:alpha, n:n, Enorm:Enorm})
    O = ones(L+1,j,C)
    R = zeros(L+1,j)
    JO = zeros1D(L+1)

    ## START SIMULATION
    ## i1 = time index. i2 = distance index
    for i1 in L:
        ## Update bulk concentrations of O and R
        for i2 in (j-2):
          O[i1+1][i2] = O[i1][i2] + DM*(O[i1][i2+1]+O[i1][i2-1]-2*O[i1][i2])
          R[i1+1][i2] = R[i1][i2] + DM*(R[i1][i2+1]+R[i1][i2-1]-2*R[i1][i2]) - km*R[i1][i2]
        JO[i1+1] = ( kf[i1+1]*O[i1+1][1] - kb[i1+1]*R[i1+1][1] ) / (1 + Dx/D*(kf[i1+1] + kb[i1+1]) )
        ## Update flux
        O[i1+1][0] = O[i1+1][1] - JO[i1+1]*(Dx/D)
        R[i1+1][0] = R[i1+1][1] + JO[i1+1]*(Dx/D) - km*R[i1+1][1]
    Z = math.eval('-n*F.*JO.*1000',{n:n, F:F, JO:JO})
    ## PLOT
    plt.plot(eta,Z)
    plt.xlabel('Overpotential (V)')
    plt.ylabel('Current density (mA/cm^2)')


def zeros1D(rows):
  ## makes 1D array of 0s
  arr = []
  for i1 in rows:
    arr[i1] = 0
  return arr

def zeros(rows,cols):
  ## makes a rowsxcol array of 0s
  arr = []
  for i1 in rows:
    for i2 in cols:
      temp[i1][i2] = 0
  return arr

def ones(rows,cols,C):
  ## Makes a rowsxcol array of C's
  arr = []
  for i1 in rows:
    for i2 in cols:
      arr[i1][i2] = C;
  return arr
