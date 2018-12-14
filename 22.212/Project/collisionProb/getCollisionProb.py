import time
from neutron import *
import matplotlib.pyplot as plt

np.random.seed(42)


class Neutron():
    def __init__(self,sideLen,reg):
        l, r = reg.x0 - reg.r, reg.x0 + reg.r
        d, u = reg.y0 - reg.r, reg.y0 + reg.r
        inCircle = False
        while not inCircle:
           self.r = np.array([l+rand()*(r-l),d+rand()*(u-d)])
           inCircle = (self.r[0]-reg.x0)**2 + (self.r[1]-reg.y0)**2 < reg.r**2

        self.mu = np.cos((2.0*rand()-1.0)*pi*0.5)
        theta = rand()*2*pi
        self.u = np.array([np.cos(theta), np.sin(theta)])


def getCollisionProb(pitch,radius,plot,numParticles,hole,fSigT_hi, \

    fSigT_lo,mSigT,verbose,startNeutronsFrom):
    sideLen  = pitch*3
    regions,surfaces = makeGeometry(pitch,radius,hole)

    start = time.perf_counter()
    pinHits, modHits = np.array([0]*9), 0
    circles = surfaces[0]
    startRegion = circles[startNeutronsFrom]

    if plot:
        ax = plt.gca(); ax.cla()
        ax.set_xlim((0,sideLen))
        ax.set_ylim((0,sideLen))

        for c in circles:
            cPlot = plt.Circle((c.x0,c.y0),c.r,color='r')
            ax.add_artist(cPlot)
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


    counter = 0
    l, r = startRegion.x0 - startRegion.r, startRegion.x0 + startRegion.r
    d, u = startRegion.y0 - startRegion.r, startRegion.y0 + startRegion.r
    for i in range(numParticles):
        inCircle = False
        while not inCircle:
           x, y = l+rand()*(r-l), d+rand()*(u-d)
           inCircle = (x-startRegion.x0)**2 + (y-startRegion.y0)**2 < startRegion.r**2

        mu = np.cos((2.0*rand()-1.0)*pi*0.5)
        theta = rand()*2*pi
        ux, uy = np.cos(theta), np.sin(theta)

        if plot: ax.plot(x,y,'bo')

        collided = False
        while not collided:
            if plot: old_x, old_y = x, y
            distTraveled,regionID,x,y,ux,uy = advance(x,y,ux,uy,mu,surfaces,regions)
            if plot:
                ax.plot(x,y,color=colors[counter%len(colors)],marker='o')
                ax.plot([old_x,x],[old_y,y],colors[counter%len(colors)]) 
                counter += 1

            if regions[regionID].name == 'mod':
                collided = weCollide(distTraveled,mSigT)
                if collided: modHits += 1

            if regions[regionID].name == 'fuel (high enr.)':
                collided = weCollide(distTraveled,fSigT_hi)
                if collided: pinHits[regionID] += 1

            if regions[regionID].name == 'fuel (low enr.)':
                collided = weCollide(distTraveled,fSigT_lo)
                if collided: pinHits[regionID] += 1



    assert(sum(pinHits)+modHits == numParticles)

    if plot: plt.show()
    
    if verbose: 
        end = time.perf_counter()
        elapsed_time = end - start
        print('Elapsed time:', '%.4f'%elapsed_time)
        print(pinHits/numParticles*100)
    return pinHits/numParticles


if __name__ == "__main__":
    plot          = False
    pitch         = 1.26
    radius        = 0.39218
    numParticles  = 10000
    hole          = True       # Pin Labels
    startNeutronsFrom = 0      #  6  7  8
    fSigT_lo = 0.27413873      #  3  4  5 
    fSigT_hi = 0.88674345      #  0  1  2 
    mSigT    = 0.26955675      
    verbose = True
    collisionProb = getCollisionProb(pitch,radius,plot,numParticles,hole,fSigT_hi,fSigT_lo,mSigT,verbose,startNeutronsFrom)

