import matplotlib
matplotlib.use('Agg')

def RHS(w, s, p):
    """
    Defines the differential equations for the system 

    Arguments:
        x :  vector of the state variables:
                  x = [t,x=r',r,phi]
        s :  affine parameter
        p :  vector of the parameters:
                  p = [rs,a,b]
                  a = p_phi, p_t = -a/b
    """
    t, x, r, phi = w
    rs, a, b = p

    # Create f = (t',x',r',phi'):
    f = [-a * 1/(1-rs/r),
         1/r**2 * (rs/2.0 + b**2/r - 3./2.0 * b**2 * rs / r**2),
         x,
         b/r**2]
    return f

def writeToFile(proper_time,solutions):
    with open('output.dat', 'w') as f:
        # Print & save the solution.
        for time, solution in zip(proper_time, solutions):
            f.write(str(time)+' '.join([" %s" % i for i in solution])+'\n')    

def solveGeodesic(RightHandSide,s,p,x0):
    
    solutions = [[0,0,0,0] for i in range(numpoints)]
    solutions[0] = x0
    rs, a, b = p

    # ODE solver parameters
    abserr = 1.0e-6
    relerr = 1.0e-6

    # loop each time step, kill loop if r is too close to rs
    for idx in range(1,int(numpoints)):
        if x0[2]>1.05*rs and x0[2]<5*r1:
            solutions[idx] = odeint(RightHandSide, x0, [s[idx-1],s[idx]], args=(p,), atol=abserr, rtol=relerr)[1]
            x0 = solutions[idx]
            #if idx>1 and idx%50 == 0:
                #print(time.time()-tstar)
                #writeToFile(s[:idx],solutions[:idx])
                #drawBitmap(rs,r1)
                #trim('output')
                #drawOnScreen()
        else:
            writeToFile(s[:idx],solutions[:idx])
            drawBitmap(rs,r1)
            trim('output')
            drawOnScreen()
            break
    
    return s[:idx], solutions[:idx]

def drawBitmap(rs,r1):
    s, t, x, r, phi = loadtxt('output.dat', unpack=True)
    figure(1, figsize=(6, 2.7))

    ax = subplot(111, polar=True)
    ax.grid(False)
    ax.axis('off')
    ax.plot(phi,r,linewidth=0.6, color='black')
    ax.set_ylim(0, r1)
    circle = Circle((0, 0), rs, transform=ax.transData._b, color='black', fill=True)
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
        width  = im.size[0]
        height = im.size[1]
        ratio = 200./90.
        if width/height > ratio:
            w1 = 0.5*width
        else:
            w1 = 0.2*width
        h1 = height - (width - w1)/ratio
        area = (w1, h1, width, height)
        cropped_im = im.crop(area)
        cropped_im.save(file+'-trimmed.png')

def drawOnScreen():
    image = PapirusImage()
    image.write('output-trimmed.png')

def run(RHS, s, p, x0):
    s, xsol=solveGeodesic(RHS,s,p,x0)

# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
from numpy import sqrt, zeros, pi, linspace, append, loadtxt, pi, abs
from pylab import figure, polar, plot, xlabel, ylabel, grid, hold, legend, title, savefig, show, subplot
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from PIL import Image, ImageChops
from papirus import PapirusImage
import time
tstart = time.time()

# Parameter values. these provide the initial velocities via line 11
rs = 1
a = 1
b = -1.755

# Initial conditions
s1 = 0
r1 = 8  # in multiples of rs
x1 = -1 # this is initial r velocity if you'd like
phi1 = 0 # this is where "around" the BH the ray starts, 0 for east, pi/2 for north, etc


sf = 10**4
numpoints = sf*20+1
s = linspace(s1,sf,numpoints)
# Pack up the parameters and initial conditions:
p = [rs, a, b]
x0 = [s1, x1, r1, phi1]

# loop over b's
for i in linspace(-3,1,31):
    run(RHS, s, [rs, a, i], x0)
trim('output')
