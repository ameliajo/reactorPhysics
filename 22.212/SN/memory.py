import matplotlib.pyplot as plt


dx = [1e-5,1e-4,1e-3,1e-2,1e-1,1,5,25]
nangles = [2e2,2e3,2e4,2e5,2e6,2e7,2e8]

for nangle in nangles:
    vec = []
    for x in dx:
        vec.append(50/x*nangle)
    plt.plot(dx,vec,label="nangle "+str('%.2E'%nangle))

plt.title('relative memory usage assuming 50 cm slab')
plt.yscale('log')
plt.ylabel('Proportional to memory usage')
plt.xlabel('dx')
plt.legend(loc='best')
plt.show()


