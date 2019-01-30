import sys
import numpy as np
import matplotlib.pyplot as plt


def xzproj_hist2d(model, bins=None, colorbar=True):
    """
    Plot a 2D histogram with the AGAMA model data
    """
    if bins is None:
        bins = [np.linspace(-20, 20, 50), np.linspace(-20, 20, 50)]
    plt.hist2d(model[:, 0], model[:, 2], bins=bins)
    if colorbar:
        plt.colorbar()
    plt.gca().set_aspect('equal')


def cyl_hist(model, bins=None, ymin=0, ymax=1000):
    """
    Plot a 1D cylindrical histogram of the inner regions of the AGAMA model data
    """
    if bins is None:
        bins = np.linspace(0, 20, 21)
    rho = (model[:, 0]**2 + model[:, 1]**2)**.5
    cyl = model[rho < 1.]
    z_cyl = np.abs(cyl[:, 2])
    plt.hist(z_cyl, bins=np.linspace(0, 20, 21))
    plt.gca().set_ylim(ymin=ymin, ymax=ymax)


def zrho_hist2d(model, bins=None, xymin=-40, xymax=40):
    """
    Plot a 2D histogram of z agains rho of the AGAMA model data
    (leaves central region empty because of automatic rounding during binning)
    """
    plt.hist2d(rho, z, bins=[np.linspace(0, 20, 50), np.linspace(0, 20, 50)])
    plt.gca().set_xlim(xmin=-xymin, xmax=xymax)
    plt.gca().set_ylim(ymin=-xymin, ymax=xymax)


def enc_mass(model):
    """
    Calculate enclosed mass plotting data from AGAMA model data
    """
    # spherical radius in kpc
    r = np.sqrt(model[:, 0]**2 + model[:, 1]**2 + model[:, 2]**2)
    R = np.logspace(-3, 3, 50)
    hist = np.histogram(r, bins=R)
    mass_prof = hist[0]*model[0, -1]
    radii = .5*(hist[1][:-1] + hist[1][1:])
    return np.cumsum(mass_prof), radii


if __name__ == "__main__":
    
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
    #cyl_hist(model)

    # plot the x-z 2d histogram
    #xzproj_hist2d(model)

    # plot the z-rho 2d histogram
    #zrho_hist2d(model)

    enc_mass(model)

    # save and/or show figure
    if len(sys.argv) > 2:
        savename = sys.argv[2]
        plt.savefig(savename)
    plt.show()
