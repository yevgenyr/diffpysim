import os
import sys
import numpy as np
import matplotlib.pyplot as plt  # Plotting
import Dans_Diffraction as dif
import scipy.special as sps
import diffpy.structure

# cf = os.path.dirname(__file__)
# sys.path.insert(0, os.path.join(cf, '..'))
dirpath = '../example_data/ran_stretchnmf'
paths = list([os.path.join(dirpath, i) for i in ["test.cif", "ZnSe_216_F-43m_c.cif", "BaTiO3_221_Pm-3m_c.cif"]])


###  https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
def G(x, alpha):
    """ Return Gaussian line shape at x with HWHM alpha """
    return np.sqrt(np.log(2) / np.pi) / alpha * np.exp(-(x / alpha) ** 2 * np.log(2))


def L(x, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / (x ** 2 + gamma ** 2)


def V(x, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(sps.wofz((x + 1j * gamma) / sigma / np.sqrt(2))) / sigma \
           / np.sqrt(2 * np.pi)


###  https://linuxtut.com/en/def7848f9194151e540c/
def lorentzian(theta, gamma):
    lorentz = 1 / (1 + (theta / gamma) ** 2)
    return lorentz


def gaussian(theta, lamb):
    sigma = 1 / 20
    gauss = (1 / (2 * np.pi * sigma)) * np.exp(-(theta ** 2 + lamb ** 2) / (2 * sigma))
    return gauss


def voigt(theta, Hg, Hl):
    z = (theta + 1j * Hl) / (Hg * np.sqrt(2.0))
    w = sps.wofz(z)
    v = (w.real) / (Hg * np.sqrt(2.0 * np.pi))
    return v


f = paths[2]

# Dans_Diffraction Examples
# Generate powder spectrum from a cif file
xtl = dif.Crystal(f)

# energy_kev = dif.fc.wave2energy(1.5498)  # 8 keV
# max_twotheta = 180

xtl.Scatter.setup_scatter('xray')  # 'xray','neutron','xray magnetic','neutron magnetic','xray resonant'
# max_wavevector = dif.fc.calqmag(max_twotheta, energy_kev)
q, I = xtl.Scatter.generate_powder(q_max=30, peak_width=0, background=0, powder_average=True)
# convolve peaks
t = np.arange(-1000, 1000, 0.1)
fwhmg = 1.5
fwhml = 1.5
v = sps.voigt_profile(t, fwhmg, fwhml)
I = np.convolve(I, v, 'same')

# add background
bkg = 0
# background = 0.01
# bkg = np.random.normal(background, np.sqrt(background), [len(I)])
# bkg = np.sin(np.ones(len(I))) * 100
I = I + bkg

plt.plot(q, I)

# convert wavevector, q=2pi/d to two-theta:
# twotheta = dif.fc.cal2theta(q, energy_kev)

# save data as csv file
# head = '%s\nx-ray powder diffraction energy=%6.3f A-1\n Q, intensity' % (xtl.name, q.max())
# np.savetxt('powder.csv', (twotheta, intensity), delimiter=',', header=head)

# # plot the spectra
# plt.figure()
# plt.plot(q, I, '-', lw=2)
# dif.fp.labels('x-ray powder diffraction E=%6.3f A-1\n%s' % (q.max(), xtl.name), 'Q', 'intensity')
# plt.show()
