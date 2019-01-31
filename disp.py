import sys
import numpy as np
import matplotlib.pyplot as pl

fname = sys.argv[1]

fil = open(fname,'r')
model = np.loadtxt(fil)

print(model.shape)

r2 = model[:,0]**2 + model[:,1]**2 + model[:,2]**2
rad = r2**.5

r = []
sv = []

rbin = np.logspace(-3,4,40)
for i in range(len(rbin)-1):
    b = (rad > rbin[i]) * (rad < rbin[i+1])
    if len(rad[b]) < 10:
        continue
    rb = np.mean(rad[b])
    svx = np.std(model[b,3])
    svy = np.std(model[b,4])
    svz = np.std(model[b,5])
    svb = (svx**2 + svy**2 + svz**2)**.5
    r.append(rb)
    sv.append(svb)
    print(rb,svb)

pl.semilogx(r,sv,'.')

pl.show()
