import h5py
import matplotlib.pyplot as plt
from plot_property import *
import numpy as np

def plot_readIn(ax1, ax2, folder, mass=1.3, flavor='charm', gridT=61, gridE=71, lst='-'):
    qhat_Qq2Qq = h5py.File('{}/qhat_Qq2Qq.hdf5'.format(folder), 'r')['Qhat-tab']
    qhat_Qg2Qg = h5py.File('{}/qhat_Qg2Qg.hdf5'.format(folder), 'r')['Qhat-tab']
    rate_Qq2Qq = h5py.File('{}/rQq2Qq.hdf5'.format(folder), 'r')['Rates-tab']
    rate_Qg2Qg = h5py.File('{}/rQg2Qg.hdf5'.format(folder), 'r')['Rates-tab']


    GeV_to_fmInv = 5.068
    fmInv_to_GeV = 0.1973

    E = np.linspace(mass*1.01, 140., gridE)
    temp = np.linspace(0.15, 0.75, gridT)
    p = np.sqrt(E**2 - mass**2)

    dpzdt = qhat_Qq2Qq[1,:,:] + qhat_Qg2Qg[1,:,:]
    kperp = qhat_Qq2Qq[2,:,:] + qhat_Qg2Qg[2,:,:]
    kpara = qhat_Qq2Qq[3,:,:] - np.square(qhat_Qq2Qq[1,:,:]) / rate_Qq2Qq  \
          + qhat_Qg2Qg[3,:,:] - np.square(qhat_Qg2Qg[1,:,:]) / rate_Qg2Qg

    if gridT == 61:
        idxT = [0, 10, 20, 30, 40]
    elif gridT == 31:
        idxT = [0, 5, 10, 15, 20]

    if gridE == 71:
        idxP = [0, 14, 28, 42, 56]
    elif gridE == 51:
        idxP = [0, 10, 20, 30, 40]

    for i in idxT:
        color = np.random.rand(3)
        ax1.plot(p, 2*kperp[:,i] * GeV_to_fmInv, lst, lw=2, label='T={} GeV'.format(temp[i]))
        #ax1.plot(p, 2*kpara[:,i] * GeV_to_fmInv, '--', lw=2, color=color)

    for j in idxP:
        color = np.random.rand(3)
        ax2.plot(temp, 2*kperp[j,:] * GeV_to_fmInv, lst, lw=2, label = 'p={} GeV/c'.format(int(p[j])))


    #ax1.set_xlim(0, 50)
    ax1.set_ylim(0, 10)
    ax1.set_xlabel(r'$p$ [GeV/c]')
    ax1.set_ylabel(r'$\hat{q}$ [GeV$^2$/fm]')
    ax1.legend(loc='upper left')
    
    ax2.set_xlim(0.15, 0.75)
    ax2.set_ylim(0, 10)
    ax2.set_xlabel(r'T [GeV/c]')
    ax2.legend(loc='upper left')
    finish()

    res = []
    tempM, EM = np.meshgrid(temp, E)
    qhat_over_T3 = 2 * kperp / temp**3
    for i in range(len(kperp)):
        for j in range(len(kperp[i])):
            dum = np.array([tempM[i,j], EM[i,j], qhat_over_T3[i,j]])
            res.append(dum)

    res = np.array(res)
    np.savetxt('gamma-table_{}.dat'.format(flavor), res, header = 'temp energy qhat_over_T3', fmt='%10.6f')


@plotfn
def Plot_DC_charm_bottom():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    plot_readIn(ax1, ax2, 'charm_build/src', mass=1.3, flavor='charm', lst='-')
    plot_readIn(ax1, ax2, 'bottom_build/src', mass=4.2, flavor='bottom', lst='--')

@plotfn
def Plot_DC_charm():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    plot_readIn(ax1, ax2, 'charm_build/src', mass=1.3, flavor='charm', lst='-')
    plot_readIn(ax1, ax2, 'charm_build_sparser/src', mass=1.3, flavor='charm', lst='--', gridT=31, gridE=51)

@plotfn
def Plot_DC_bottom():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    plot_readIn(ax1, ax2, 'bottom_build/src', mass=4.2, flavor='bottom', lst='-')
    plot_readIn(ax1, ax2, 'bottom_build_sparser/src', mass=4.2, flavor='bottom', lst='--', gridT=31, gridE=51)



#Plot_DC_charm_bottom()
#Plot_DC_charm()
Plot_DC_bottom()

