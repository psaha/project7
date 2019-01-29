import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    print("Input filename")
    sys.exit(1)
print(filename)

with open(filename, 'r') as f:
    model = np.loadtxt(f)
    
# select slices
rho = (model[:, 0]**2 + model[:, 1]**2)**.5
cyl = model[rho < 1.]
z_cyl = np.abs(cyl[:, 2])
z = np.abs(model[:, 2])

# print some numbers
print("Particle masses: {:e} Msol".format(model[0, -1]))
print("Total mass: {:e} Msol".format(np.sum(model[:, -1])))
print("Selected cylinder: particles within 1 kpc in xy plane")
print("Total mass in cylinder: {:e} Msol".format(np.sum(cyl[:, -1])))

# plot the cylindrical histogram
#plt.hist(z_cyl, bins=np.linspace(0, 20, 21))
#plt.gca().set_ylim(ymin=0, ymax=1000)

# plot the x-z 2d histogram
plt.hist2d(model[:, 0], model[:, 2], bins=[np.linspace(-20, 20, 50), np.linspace(-20, 20, 50)])
plt.colorbar()
plt.gca().set_aspect('equal')

# plot the z-rho 2d histogram (leaves central region empty because of automatic rounding during binning)
#plt.hist2d(rho, z, bins=[np.linspace(0, 20, 50), np.linspace(0, 20, 50)])
#plt.gca().set_xlim(xmin=-40, xmax=40)
#plt.gca().set_ylim(ymin=-40, ymax=40)


# save and/or show figure
if len(sys.argv) > 2:
    savename = sys.argv[2]
    plt.savefig(savename)
plt.show()
