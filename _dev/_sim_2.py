"""Simulate PDF of expanding lattices in multiple phases."""
from diffpysim import module_path

import os

import numpy as np
import matplotlib.pyplot as plt

from pyobjcryst import loadCrystal
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe


class Plot():
    def __init__(self, array):
        self.array = array
        self.fig = None
        self.ax = None

    def plot(self, x, y):
        if not self.fig:
            self.fig, self.ax = plt.subplots(1, 1, num="simulated results")
        self.ax.plot(x, y)
        self.ax.set_xlabel(r"$r (\AA)$")
        self.ax.set_ylabel(r"$G (\AA^{-2})$")
        self.ax.legend(self.labels)

    def set_lables(self, vals):
        self.labels = [r"a,b,c = {} $\AA$".format(a.round(2)) for a in np.squeeze(vals)]

    def show(self):
        if self.fig:
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()
        else:
            plt.show()

class Sim():
    def __init__(self):
        self.experiment = 'ran_stretchnmf'
        self.datpath = os.path.join(module_path, '..', 'example_data')
        self.cifdirpath = os.path.join(self.datpath, self.experiment)
        self.files = os.listdir(self.cifdirpath)

        # =================================================================
        self.datname = os.path.join(self.cifdirpath, "test.gr")
        self.compounds = [('ZnSe', "ZnSe_216_F-43m_c.cif"),
                          ("BaTiO3", "BaTiO3_221_Pm-3m_c.cif")]
        for comp in self.compounds:  # add attributes as compounds
            setattr(self, comp[0], os.path.join(self.cifdirpath, comp[1]))

        self.range = (5.62, 5.82, 0.1)
        self.latice_names = ('a', 'b', 'c')
        self.qdamp = 0.03
        self.qbroad = 0.02
        self.qmin = 0.1
        self.qmax = 10
        self.rmin = 0
        self.rmax = 30
        self.rstep = 0.01
        self.xvector = np.arange(self.rmin, self.rmax, self.rstep)

        self.ratio = 0.5
        # =================================================================
        self.gens = list(None for i in self.compounds)
        self.phases = list(None for i in self.compounds)
        self.gen_names = list(comp[0] for comp in self.compounds)

    def makeRecipe(self):
        """Create a fitting recipe for crystalline PDF data."""

        # The Profile

        profile = Profile()
        parser = PDFParser()
        parser.parseFile(self.datname)
        profile.loadParsedData(parser)
        profile.setCalculationRange(xmax=self.rmax)

        # The ProfileGenerator
        for i, comp in enumerate(self.compounds):
            self.gens[i] = PDFGenerator(self.gen_names[i])
            self.gens[i].setStructure(loadCrystal(eval(f"self.{comp[0]}")))

        # The FitContribution
        contribution = FitContribution("mixture")
        for i in range(len(self.gen_names)):
            contribution.addProfileGenerator(self.gens[i])
        contribution.setProfile(profile, xname="r")

        # Write the fitting equation.
        contribution.setEquation(f"scale * {'+'.join(self.gen_names)}")

        # Make the FitRecipe and add the FitContribution.
        recipe = FitRecipe()
        recipe.addContribution(contribution)

        ## Configure the fit variables
        # Start by configuring the scale factor and resolution factors.
        # We want the sum of the phase scale factors to be 1.
        recipe.newVar("scale", self.ratio)
        for i in range(len(self.gen_names)):  # TODO - adjust for more than a two variable
            if i == 0:
                recipe.constrain(self.gens[i].scale, "scale")
            else:
                recipe.constrain(self.gens[i].scale, f"{i} - scale")

        # We also want the resolution factor to be the same on each.
        for var in ['qdamp', 'qbroad']:
            recipe.newVar(var, eval(f"self.{var}"))
            for i in range(len(self.gen_names)):
                setattr(self.gens[i], var, 0)  # setup
                recipe.constrain(eval(f"self.gens[{i}].{var}"), var)

        # Vary the gloabal scale as well.
        # recipe.addVar(contribution.scale, 1)

        # First phase parameters
        for i, comp in enumerate(self.compounds):
            self.phases[i] = self.gens[i].phase
            for par in self.phases[i].sgpars:
                recipe.addVar(par, name=par.name + f"_{comp[0]}")
            recipe.addVar(self.gens[i].delta2, name=f"delta2_{comp[0]}")

        # Give the recipe away so it can be used!
        self.recipe = recipe
        return self.recipe

    def plotResults(self):
        """Plot the results contained within a refined FitRecipe."""

        r = self.recipe.mixture.profile.x
        g = self.recipe.mixture.profile.y
        gcalc = self.recipe.mixture.profile.ycalc
        # diffzero = -0.8 * max(g) * np.ones_like(g)
        # diff = g - gcalc + diffzero

        import pylab
        pylab.plot(r, g, 'bo', label="G(r) Data")
        # pylab.plot(r, gcalc, 'r-', label="G(r) Fit")
        # pylab.plot(r, diff, 'g-', label="G(r) diff")
        # pylab.plot(r, diffzero, 'k-')
        pylab.xlabel(r"$r (\AA)$")
        pylab.ylabel(r"$G (\AA^{-2})$")
        pylab.legend(loc=1)

        pylab.show()

    def main(self):
        self.recipe = self.makeRecipe()
        self.plotResults()
        return self.recipe

if __name__ == '__main__':
    # gen = Sim().main()
    gen = Sim().main()
