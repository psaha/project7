import sys
from agama_reader import enc_mass


def cosmo_cfunc(M200,h):
    #From Dutton & Maccio 2014:
    c = 10.**(0.905 - 0.101 * (np.log10(M200*h)-12.))
    return c

def corenfw_mass(r,M200,c,oden,tSF,Rhalf,rhocrit,eta,kappa):
    #Assumes input arrays in Msun, kpc:
    G = 6.67e-11
    kpc = 3.086e19
    Msun = 1.989e30
    Gyr = 365.*24.*60.*60.*1e9

    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    print 'coreNFW scale length:', rs

    Mrs = M200 * gcon * (np.log(2.0)-0.5)
    tdyn_rs = \
        2.*np.pi*np.sqrt((rs*kpc)**3./G/(Mrs*Msun))/Gyr

    rc = eta * Rhalf

    x = r/rc
    f = (np.exp(x) - np.exp(-x))/(np.exp(x)+np.exp(-x))
    xx = tSF/tdyn_rs * kappa
    if (tSF > 0.0):
        n = (np.exp(xx) - np.exp(-xx))/(np.exp(xx)+np.exp(-xx))
    else:
        n = -tSF
    my_manal = manal*f**n

    my_rhoanal = rhoanal*f**n + \
        1.0/(4.*np.pi*r**2.*rc)*manal*(1.0-f**2.)*n*f**(n-1.0)

    return my_manal


#Simple program to make a plot of the MW dark matter cumulative mass profile.
if __name__ == "__main__":
    #Import plots library:
    import numpy as np
    import pylab as plt
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    #Parameters / options:
    include_anal_abund = 'no'
    include_anal_range = 'no'
    use_fillbet = 'no'
    include_core = 'yes'

    # Read in data:
    f = open('Wegg_DM_mass_profile_MW.txt','r')
    data_W = np.genfromtxt(f)
    f.close()
    f = open('Portail17.txt','r')
    data_P = np.genfromtxt(f)
    f.close()
    f = open('ColeBinney.txt','r')
    data_CB = np.genfromtxt(f)
    f.close()
    f = open('Kafle14.txt','r')
    data_K = np.genfromtxt(f)
    f.close()

    #Mstar-Mhalo relation:
    f = open('mstar-mhalo-field.txt','r')
    data = np.genfromtxt(f)
    f.close()

    #Work out abundance matching MW mass and DM profile:
    M200lowarr = data[:,0][np.where(data[:,1])]
    Mstarlow = data[:,1][np.where(data[:,1])]
    M200higharr = data[:,0][np.where(data[:,2])]
    Mstarhigh = data[:,2][np.where(data[:,2])]
    MstarMW_low = (1.84-0.7)*1e10 + (4.6-0.3)*1e10
    MstarMW_high = (1.84+0.7)*1e10 + (4.6+0.3)*1e10
    M200MW_high = np.interp(MstarMW_high,Mstarlow,M200lowarr) 
    M200MW_low = np.interp(MstarMW_low,Mstarlow,M200lowarr)
    print 'Stellar mass MW [1e10 Msun]:',\
        MstarMW_low/1e10, MstarMW_high/1e10
    print 'Mass MW (abund match):',\
        M200MW_low, M200MW_high
    h = 0.7
    clow = cosmo_cfunc(M200MW_low,h)
    clow = 10.0**(np.log10(clow)-0.1)
    chigh = cosmo_cfunc(M200MW_high,h)
    chigh = 10.0**(np.log10(chigh)+0.1)
    print 'concentration par range:',clow, chigh

    #Analytic profiles:
    ranal = np.logspace(-1,3,500)
    oden = 200.0
    tSF_cusp = 0.01
    Rhalf_cusp = 1.0
    tSF_core = -1
    Rhalf_core = 1.0
    rhocrit = 135.05
    eta = 1.75
    kappa = 0.04
    manal_low_cusp = corenfw_mass(ranal,M200MW_low,clow,\
                                  oden,tSF_cusp,Rhalf_cusp,rhocrit,eta,kappa)
    manal_high_cusp = corenfw_mass(ranal,M200MW_high,chigh,\
                                   oden,tSF_cusp,Rhalf_cusp,rhocrit,eta,kappa)
    manal_low_core = corenfw_mass(ranal,M200MW_low,clow,\
                                  oden,tSF_core,Rhalf_core,rhocrit,eta,kappa)
    manal_high_core = corenfw_mass(ranal,M200MW_high,chigh,\
                                   oden,tSF_core,Rhalf_core,rhocrit,eta,kappa)

    #Range of halo masses to consider:
    manal_range = [5e11,1e12,5e12]

    #Additional data:
    Posti_Helmi = 1.37e11
    Posti_Helmi_err = 0.12e11
    Posti_Helmi_r = 20.0

    Bovy_streams = 1.1e11
    Bovy_streams_err = 0.1e11
    Bovy_streams_r = 20.0

    Watkins_r = [21.1,39.5]
    Watkins = [0.22e12,0.44e12]
    Watkins_lerr = [0.03e12,0.06e12]
    Watkins_herr = [0.04e12,0.07e12]

    #Portail data is *density* in Msun pc^-3; need to convert
    #to M(<r):
    M_P = np.zeros(len(data_P))
    M_P_low = np.zeros(len(data_P))
    M_P_high = np.zeros(len(data_P))
    r_P = data_P[:,0]*1000.0
    r_P_use = np.zeros(len(data_P))
    den_P = data_P[:,1]
    den_P_low = data_P[:,1]+data_P[:,2]
    den_P_high = data_P[:,1]+data_P[:,3]
    for i in range(len(data_P)-1):
        M_P[i+1] = M_P[i] + \
            4.0/3.0*np.pi*(r_P[i+1]**3.0-r_P[i]**3.0)*den_P[i]
        M_P_low[i+1] = M_P_low[i] + \
            4.0/3.0*np.pi*(r_P[i+1]**3.0-r_P[i]**3.0)*den_P_low[i]
        M_P_high[i+1] = M_P_high[i] + \
            4.0/3.0*np.pi*(r_P[i+1]**3.0-r_P[i]**3.0)*den_P_high[i]
        r_P_use[i+1] = r_P[i+1]/1000.0
    M_P = M_P[1:]
    M_P_low = M_P_low[1:]
    M_P_high = M_P_high[1:]
    r_P_use = r_P_use[1:]

    #ColeBinney is vc(r). Convert to M(<r):
    Gsi = 6.67e-11
    Msun = 1.989e30
    kpc = 3.086e19
    kms = 1e3
    Guse = Gsi * Msun / kpc / kms**2.0
    M_CB = data_CB[:,1]**2.0*data_CB[:,0]/Guse
    M_CB_low = (data_CB[:,1]+data_CB[:,2])**2.0*data_CB[:,0]/Guse
    M_CB_high = (data_CB[:,1]+data_CB[:,3])**2.0*data_CB[:,0]/Guse

    #Make the plot:
    figsize = 8
    figx = figsize
    figy = figsize
    myfontsize = 30
    mylinewidth = 3
    mylinewidth2 = 4

    #Use LaTeX Fonts: 
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    #Set thick axes: 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)

    #Make nice tick marks:
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    plt.loglog()

    plt.xlabel(r'$r\,[{\rm kpc}]$',fontsize=myfontsize)
    plt.ylabel(r'$M_{\rm DM}(<r)\,[{\rm M}_\odot]$',\
                   fontsize=myfontsize)

    plt.fill_between(data_W[:,0],data_W[:,1]+data_W[:,2],\
                         data_W[:,1]+data_W[:,3],\
                         facecolor='0.65',\
                         edgecolor='none',zorder=1)
    plt.plot(data_W[:,0],data_W[:,1],'k',linewidth=mylinewidth,\
                 label=r'Wegg et al. 2018 $\mid$ RR-Lyrae',zorder=2)

    plt.fill_between(r_P_use,M_P_low,\
                         M_P_high,\
                         facecolor='r',\
                         edgecolor='none',zorder=1,alpha=0.5)
    plt.plot(r_P_use,M_P,'r',linewidth=mylinewidth,\
                 label=r'Portail et al. 2017 $\mid$ Inner Gal.',zorder=2)

    plt.fill_between(data_CB[:,0],M_CB_low,\
                         M_CB_high,\
                         facecolor='b',\
                         edgecolor='none',zorder=1,alpha=0.5)
    plt.plot(data_CB[:,0],M_CB,'b',linewidth=mylinewidth,\
                 label=r'Cole \& Binney 2017 $\mid$ Giant stars',zorder=2)


    plt.errorbar(Posti_Helmi_r,Posti_Helmi,\
                 yerr=Posti_Helmi_err,fmt='o',\
                 markersize=5,color='purple',ecolor='purple',\
                 label=r'Posti \& Helmi 2018 $\mid$ GCs',\
                 zorder=3)

    plt.errorbar(Bovy_streams_r,Bovy_streams,\
                 yerr=Bovy_streams_err,fmt='o',\
                 markersize=5,color='magenta',ecolor='magenta',\
                 label=r'Bovy et al. 2016 $\mid$ Streams',\
                 zorder=3)

    plt.errorbar(Watkins_r,Watkins,\
                 yerr=[Watkins_lerr,Watkins_herr],fmt='o',\
                 markersize=5,color='green',ecolor='green',\
                 label=r'Watkins et al. 2018 $\mid$ GCs',\
                 zorder=3)

    plt.fill_between(data_K[:,0],data_K[:,1]+data_K[:,2],\
                         data_K[:,1]+data_K[:,3],\
                         facecolor='green',alpha=0.5,\
                         edgecolor='none',zorder=1)
    plt.plot(data_K[:,0],data_K[:,1],'g',linewidth=mylinewidth,\
                 label=r'Kafle et al. 2014 $\mid$ halo stars',zorder=2)

    if (include_anal_abund == 'yes'):
        plt.fill_between(ranal,manal_low_cusp,manal_high_cusp,\
                         facecolor='orange',alpha=0.75,\
                         edgecolor='none',zorder=0)
    if (include_core == 'yes'):
        plt.fill_between(ranal,manal_low_core,manal_high_core,\
                             facecolor='orange',alpha=0.75,\
                             edgecolor='none',zorder=0)

    if (include_anal_range == 'yes'):
        for i in range(len(manal_range)):
            M200MW = manal_range[i]
            clow = cosmo_cfunc(M200MW,h)
            clow = 10.0**(np.log10(clow)-0.1)
            chigh = cosmo_cfunc(M200MW,h)
            chigh = 10.0**(np.log10(chigh)+0.1)
            manal_low_cusp = corenfw_mass(ranal,M200MW,clow,\
                                  oden,tSF_cusp,Rhalf_cusp,rhocrit,eta,kappa)
            manal_high_cusp = corenfw_mass(ranal,M200MW,chigh,\
                                   oden,tSF_cusp,Rhalf_cusp,rhocrit,eta,kappa)
            manal_low_core = corenfw_mass(ranal,M200MW,clow,\
                                  oden,tSF_core,Rhalf_core,rhocrit,eta,kappa)
            manal_high_core = corenfw_mass(ranal,M200MW,chigh,\
                                   oden,tSF_core,Rhalf_core,rhocrit,eta,kappa)

            if (use_fillbet == 'yes'):
                plt.fill_between(ranal,manal_low_cusp,manal_high_cusp,\
                         facecolor='brown',alpha=0.75,\
                         edgecolor='none',zorder=0)
                if (include_core == 'yes'):
                    plt.fill_between(ranal,manal_low_core,manal_high_core,\
                             facecolor='brown',alpha=0.75,\
                             edgecolor='none',zorder=0)
            else:
                plt.plot(ranal,(manal_low_cusp+manal_high_cusp)/2.0,\
                         linewidth=2,linestyle='--',color='black',\
                         zorder=0)

    # compare with AGAMA model data
    if len(sys.argv) > 1:
        filename = sys.argv[1:]
    else:
        print("Input filename")
        sys.exit(1)
    print(filename)

    model = []
    for i, fyle in enumerate(filename):
        with open(fyle, 'r') as f:
            model.append(np.loadtxt(f))

        Menc, radii = enc_mass(model[i])
        plt.plot(radii, Menc, lw=2, label="AGAMA{}".format(i))

    # limit between core and cusp
    limrads = np.logspace(-3, 3, 500)
    plt.plot(limrads, 1e9*limrads**3, lw=2, ls='--', color='black')
                    
    plt.ylim([1e6,3e12])
    plt.xlim([0.1,300])
    
    plt.legend(loc='lower right',fontsize=15)
    # plt.savefig('dmhalo_massr.pdf')
    plt.show()
