import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import pi
from math import exp
import matplotlib.gridspec as gridspec

def getSAB(a,b):
    return (4.0*pi*a)**-0.5 * exp(-(a+b)**2/(4.0*a))



fig = plt.figure()
ax2 = fig.add_subplot(111)
plt.xlabel("Final Energy E' (eV)")
plt.ylabel('Free Gas S(a,b)')
plt.title('Free Gas S(a,b) for room temperature, E = 1 eV')

#ax2 = fig.add_subplot(121)
ax = fig.add_axes([0.93,0.15,0.05,0.7])
#gs = gridspec.GridSpec(4,1)
#ax.set_position(gs[0:2].get_position(fig))



E = 1.0
kbT = 0.025
betas = np.linspace(-E/kbT,10,500)
alphas = [5e-1]+list(np.linspace(1.0,20,10))

allSAB = []
Eprimes = [beta*kbT+E for beta in betas]

from pylab import *
cmap = cm.get_cmap('tab10', len(alphas))    # PiYG
for a,alpha in enumerate(alphas):
    rgb = cmap(a)[:3]
    sab = [getSAB(alpha,beta) for beta in betas]
    allSAB.append(sab)
    ax2.plot(Eprimes,sab,label="alpha = "+str(alpha),color=matplotlib.colors.rgb2hex(rgb),linewidth=2)

#plt.legend(loc='best')
norm = matplotlib.colors.Normalize(vmin=alphas[0], vmax=alphas[-1])

cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical',label='alphas')
plt.savefig('temp.png',bbox_inches='tight',pad_inches=0.5)


