import matplotlib
matplotlib.use('Agg')

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

def write_to_file(proper_time,solutions):
    with open('output.dat', 'w') as f:
        # Print & save the solution.
        for time, solution in zip(proper_time, solutions):
            f.write(str(time)+' '.join([" %s" % i for i in solution])+'\n')    

def solveGeodesic(RightHandSide,tau,p,x0):
    
    solutions = [[0,0,0] for i in range(numpoints)]
    solutions[0] = x0

    # ODE solver parameters
    abserr = 1.0e-8
    relerr = 1.0e-6

    # loop each time step, kill loop if r is too close to rs
    for idx in range(1,int(numpoints)):
        if x0[1]>1.05*rs and x0[1]<10*rs:
            #print(tau[idx],x0)
            solutions[idx] = odeint(RightHandSide, x0, [tau[idx-1],tau[idx]], args=(p,), atol=abserr, rtol=relerr)[1]
            x0 = solutions[idx]
        else:
            # if we have reached the horizon, simply fill in the rest of the solution array with previous values
            solutions[idx] = solutions[idx-1]
    
    return solutions

def drawBitmap(rs):
    tau, t, r, phi = loadtxt('output.dat', unpack=True)
    figure(1, figsize=(6, 2.7))

    ax = subplot(111, polar=True)
    ax.grid(False)
    ax.axis('off')
    ax.plot(phi,r,linewidth=0.8, color='black')
    ax.set_ylim(0, 6)
    initial = Circle((r[0], phi[0]), 0.02, transform=ax.transData._b, color='red', fill=True)
    circle = Circle((0, 0), rs, transform=ax.transData._b, color='black', fill=True)
    ax.add_artist(initial)
    ax.add_artist(circle)

    savefig('output.png', dpi=800)


def trim(file):
    im = Image.open(file+'.png')
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        im = im.crop(bbox)
        im.save(file+'-trimmed.png')

def drawOnScreen():
    image = PapirusImage()
    image.write('output-trimmed.png')

# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
from numpy import sqrt, zeros, pi, linspace, append, loadtxt, pi
from pylab import figure, polar, plot, xlabel, ylabel, grid, hold, legend, title, savefig, show, subplot
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from PIL import Image, ImageChops
from papirus import PapirusImage

# Parameter values. these provide the initial velocities via line 11
rs = 1
a = -2          # changing sign seems to rotate p_phi by pi. value doesn't seem to do much
b = 2.4       # sign changes direction from the BH but still same initial p_phi it seems
                # higher abs value increases deviation from straight to BH
                # b = 2.598 seems to be close to maximum
                # very small values just send it straight in
                # b [eps,2.598]



# Initial conditions
t1 = 0
r1 = 8  # in multiples of rs
phi1 = 0 # this is where "around" the BH the ray starts, 0 for east, pi/2 for north, etc


tf = 10**4
numpoints = tf*20+1
tau = linspace(t1,tf,numpoints)
# Pack up the parameters and initial conditions:
p = [rs, a, b]
x0 = [t1, r1, phi1]

xsol=solveGeodesic(RHS,tau,p,x0)
write_to_file(tau,xsol)
drawBitmap(rs)
trim('output')
drawOnScreen()
