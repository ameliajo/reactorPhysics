import matplotlib.pyplot as plt 
import numpy as np
from scipy import interpolate


with open('q3_data.txt') as f:
    lines = f.readlines()
energy = [0.0]*len(lines)
xs     = [0.0]*len(lines)
for i,line in enumerate(lines):
    energy[i] = float(line.split(',')[0])
    xs[i]     = float(line.split(',')[1])



RI_integrand = [xs[i]/energy[i] for i in range(len(energy))]


ranges = [[0.50,1.e4],\
          [0.01,0.10],\
          [0.10,1.00],\
          [1.00,6.00],\
          [6.00,10.0],\
          [10.0,25.0],\
          [25.0,50.0],\
          [50.0,100.]]




# This is if you just walk through and pick out the points that 
# fit in the range
for range_ in ranges:
    xList = []
    yList = []
    for i in range(len(energy)):
        if (energy[i] > range_[0] and energy[i] <= range_[1]):
            xList.append(energy[i])
            yList.append(RI_integrand[i])

    print(np.trapz(yList,xList))
print('')


# This is if you interpolate 
f = interpolate.interp1d(energy, RI_integrand)

RI = []

for range_ in ranges:
    xlist = np.linspace(range_[0],range_[1],5000)
    ylist = list(f(xlist))
    RI.append(np.trapz(ylist,x=xlist))
    plt.fill_between([range_[0],range_[0],\
                      range_[1],range_[1]],\
                      [0,RI[-1],RI[-1],0],alpha=0.7)
    print(range_,RI[-1])

plt.plot(energy,xs)
minVal = min([x[0] for x in ranges])
maxVal = max([x[1] for x in ranges])
plt.xlim([minVal,maxVal])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('Cross Section (barn)')
plt.title('Q3: U-238 Pointwise Radiative Capture and Resonance Integrals')
plt.show()








