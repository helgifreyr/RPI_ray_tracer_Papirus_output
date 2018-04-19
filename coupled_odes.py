def RHS(x, tau, p):
    """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        x :  vector of the state variables:
                  x = [t,r,phi]
        tau :  proper time
        p :  vector of the parameters:
                  p = [rs,a,b]
                  a = p_phi, p_t = -a/b
    """
    t, r, phi = x
    rs, a, b = p

    # Create f = (t',r',phi'):
    f = [-a/b * 1/(1-rs/r),
         (a/b) * sqrt(1 - b**2 * (1-rs/r)/r**2),
         -a/r**2]
    return f

# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
from numpy import sqrt, zeros, pi, linspace, append

# Parameter values. these provide the initial velocities via line 11
rs = 1
a = -2          # changing sign seems to rotate p_phi by pi. value doesn't seem to do much
b = 2.598       # sign changes direction from the BH but still same initial p_phi it seems
                # higher abs value increases deviation from straight to BH
                # b = 2.598 seems to be close to maximum
                # very small values just send it straight in
                # b [eps,2.598]



# Initial conditions
t1 = 0
r1 = 8  # in multiples of rs
phi1 = 0 # this is where "around" the BH the ray starts, 0 for east, pi/2 for north, etc

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
tf = 10**4
numpoints = tf*20+1

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
tau = linspace(t1,tf,numpoints)

# Pack up the parameters and initial conditions:
p = [rs, a, b]
x0 = [t1, r1, phi1]

xsol = [[0,0,0] for i in range(numpoints)]
xsol[0] = x0

# loop each time step, kill loop if r is too close to rs
for idx in range(1,int(numpoints)):
    if x0[1]>1.05*rs and x0[1]<10*rs:
        #print(tau[idx],x0)
        xsol[idx] = odeint(RHS, x0, [tau[idx-1],tau[idx]], args=(p,), atol=abserr, rtol=relerr)[1]
        x0 = xsol[idx]
    else:
        # if we have reached the horizon, simply fill in the rest of the solution array with previous values
        xsol[idx] = xsol[idx-1]
#       print(tau[idx],xsol[idx])

    
# Call the ODE solver.
#xsol = odeint(RHS, x0, tau, args=(p,),
#             atol=abserr, rtol=relerr)

with open('output.dat', 'w') as f:
    # Print & save the solution.
    for tau1, x1 in zip(tau, xsol):
        f.write(str(tau1)+' '.join([" %s" % i for i in x1])+'\n')
