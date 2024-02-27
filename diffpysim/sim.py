"""Simulate PDF of expanding lattices in multiple phases."""

import os
import yaml
import pandas as pd
import numpy as np
from tqdm import tqdm

# xrd
from diffpysim.cifmap import CifMap
import scipy.special as sps

# pdf
from pyobjcryst.crystal import Crystal
from pyobjcryst import loadCrystal
from diffpy.srfit.pdf import PDFGenerator
from diffpy.srfit.pdf.characteristicfunctions import sphericalCF

from diffpysim.plot import Plot
from diffpysim.userconfig import UserConfig


class Sim(UserConfig):
    def __init__(self):
        super().__init__()

        if not self.DUMP:
            self.DUMP_CSV = False
            self.DUMP_TXT = False
        if self.parent_path.strip()[0] == '~':
            self.parent_path = os.path.expanduser(self.parent_path.strip())
        else:
            self.parent_path = os.path.abspath(self.parent_path)
        os.chdir(self.parent_path)
        self.cifdirpath = os.path.join(self.parent_path, self.cif_bank_dirname)
        self.uiso_step = (self.uiso_max - self.uiso) / (self.steps - 1)  # for step uiso

        # =================================================================
        # code dumpster
        # self.files = os.listdir(self.cifdirpath)  # use if all files are asked
        # components and file names
        # [comp_name, {comp_filename, weight, lattice_attributes, strain}]
        # for comp in self.compounds.items():
        #     assert comp[1]['lat_attrs'] == list(self.compounds.items())[0][1]['lat_attrs']

        # ----- define structural configs

        # UC1 -- self.uiso_max = self.uiso  ,  max_weight_shift = 0
        # UC2 -- self.uiso_max (=0.01)> self.uiso  ,  max_weight_shift = 0
        # UC3 -- self.uiso_max (=0.01) > self.uiso  ,  max_weight_shift = 0.5, -0,5

        # 'Ru': {
        #     'cifname': "Ru_mp-8639_conventional_standard.cif",
        #     'weight': 0.5,
        #     'max_weight_shift': -0.5,
        #     'lat_attrs': self.lat_attrs,
        #     'max_strain': (0.005, 0.0050, 0.0050),
        #     'size': 0.0
        # },
        #
        # 'RuO2_t': {
        #     'cifname': "RuO2_tetragonal.cif",
        #     'weight': 0.5,
        #     'max_weight_shift': -0.5,
        #     'lat_attrs': self.lat_attrs,
        #     'max_strain': (0.020, 0.020, 0.020),
        #     'size': 0.0
        # },
        #
        # 'LiRuO2': {
        #     'cifname': "LiRuO2_mp-28254_conventional_standard.cif",
        #     'weight': 0.0,
        #     'max_weight_shift': 1,
        #     'lat_attrs': self.lat_attrs,
        #     'max_strain': (0.010, 0.010, 0.010),
        #     'size': 0.0
        # },

        # =================================================================

        #  =========================== CREATE MIXED COMP. ===========================

        # ----- name
        self.mixname = 'mix'
        for name, config in self.compounds.items():
            self.mixname += "___"
            weight = str(config['weight'])
            self.mixname += str(name)
            self.mixname += str("-w_" + weight)

        # ----- set avg initial strain for mix - used as metadata only
        self.avg_max_strain = tuple(self._get_avg_strains().round(3).tolist())
        # NOTE: self.avg_max_strain does not impact the actual linear combination of the simulated data. used only as metadata
        self.mix_compounds = {
            self.mixname: {
                'lat_attrs': self.lat_attrs,
                'max_strain': self.avg_max_strain
            }
        }

        # =========================== SINK (plot, save, etc.) ===========================
        self.comps_to_sink = list()
        self.comps_to_sink.extend(
            list(self.compounds.keys()))  # adding all active compounds (can choose to be selective with names)
        self.comps_to_sink.append(self.mixname)  # adding the mix compound

        # =========================== SINK (plot, save, etc.) ===========================

        # =========================================================================

        # =========================== FINALIZE SETUPS ===========================
        # ----- add comp name attributes to class
        for comp in self.compounds.items():
            compname = comp[0]
            setattr(self, compname, os.path.join(self.cifdirpath, comp[1]['cifname']))

        # ----- extract particle
        self.psize = dict(
            [(comp, self.compounds[comp]['size']) for comp in self.compounds if comp in self.comps_to_sink])
        self.psize.update({self.mixname: list(self.psize.values())})

        # ----- setup weights
        self.weights = np.array(list(comp[1]['weight'] for comp in self.compounds.items()))
        if self.NORM_WEIGHTS:
            self._norm_weights(stretch_1_0=False)

        # ----- generate weights list from compounds
        if self.SHIFT_WEIGHTS:
            self.weight_deltas = np.array(
                list(comp[1]['max_weight_shift'] / self.steps for comp in self.compounds.items()))

        # ----- memorize to metadata
        self.sim_config = dict(
            steps=self.steps,
            qdump=self.qdamp,
            qbroad=self.qbroad,
            qmin=self.qmin,
            qmax=self.qmax,
            uiso=self.uiso,
            uiso_max=self.uiso_max,
            rmin=self.rmin,
            rmax=self.rmax,
            rstep=self.rstep,
            psize=self.psize,
        )

        # ----- setup main containers
        self.gens = list(None for i in self.compounds)
        self.phases = list(None for i in self.compounds)
        self.gen_names = list(comp[0] for comp in self.compounds)
        self.pdfdict = dict(list([(comp, dict()) for comp in self.compounds]))
        self.xrddict = dict(list([(comp, dict()) for comp in self.compounds]))
        self.combos_df = pd.DataFrame()

    @staticmethod
    def set_generator(genname: str, stru: Crystal):
        """
        sets up a generator
        Parameters
        ----------
        genname: str
            the name of the generator
        stru: Crystal
            the structure object

        Returns
        -------
        PDFgenerator with the name 'gname', having the structure 'stru'
        """
        gen = PDFGenerator(genname)
        gen.setStructure(stru)
        return gen

    def set_gen_config(self, gen: PDFGenerator):
        """ sets the generator's configs"""
        gen.qdamp = self.qdamp
        gen.qbroad = self.qbroad
        gen.setQmax(self.qmax)
        gen.setQmin(self.qmin)
        self.set_uiso(gen, self.uiso)

    def set_uiso(self, gen, uiso):
        atoms = gen.phase.getScatterers()
        biso = self._get_biso_from_uiso(uiso)
        for atom in atoms:
            atom.Biso.setValue(biso)

    @staticmethod
    def _get_biso_from_uiso(uiso):
        return 8 * np.pi ** 2 * uiso

    def get_lat_vals(self, gen: PDFGenerator, func=None, *fargs, **fkwargs):
        """
        generates list of lattice coordinates to the lat_attrs
        Parameters
        ----------
        gen: PDFGenerator
            produced generator
        func:
            a function that can be applied
        fargs:
            args that relate to the function
        fkwargs
            kwargs that relate to the function
        Returns
        -------
        list of lattice coordinates following lat_attrs
        """
        lat_val_dict = dict()
        lat = gen.phase.getLattice()
        lat_attrs = self.compounds[gen.name]['lat_attrs']
        max_strain = self.compounds[gen.name]['max_strain']
        for i, attr in enumerate(lat_attrs):
            min_val = getattr(lat, attr).value
            max_val = min_val + min_val * max_strain[i]
            lat_span = np.linspace(min_val, max_val, num=self.steps, endpoint=True)
            # TODO - test decorating map over linspace
            if func:
                lat_span = map(lambda m: func(m, *fargs, **fkwargs), lat_span)
                # map(lambda m: func(m, *args, **kwargs), filter(None, [1,2,3,4]))
            lat_val_dict.update({attr: lat_span})
        lat_df = pd.DataFrame(lat_val_dict).round(3)
        return lat_df

    def generate_nano_particle_envelop(self, compname):
        envelop = list()
        psize = self.psize[compname]
        if psize:
            for r in self.xvector:
                envelop.append(float(sphericalCF(r, psize)))
            nanop_envelop = np.array(envelop)
            return nanop_envelop
        else:
            return 1

    # ================= UTIL FUNCS =====================

    def _norm_weights(self, stretch_1_0=False):
        if max(self.weights) == min(self.weights):
            self.weights = list([1 / len(self.weights) for x in self.weights])
        elif stretch_1_0:
            self.weights = (self.weights - min(self.weights)) / ((max(self.weights) - min(self.weights)))
        else:
            self.weights = np.array([w / sum(self.weights) for w in self.weights])

    def _get_avg_strains(self):
        """ gets the average strain based on the strain of the components """
        strains = list()
        avg_strain = list()
        for name, config in self.compounds.items():
            strains.append(np.array(config['max_strain']) * config['weight'])
        strains = np.asarray(strains)
        for s in strains.T:
            avg_strain.append(np.sum(s))
        avg_strain = np.array(avg_strain)
        return avg_strain

    @staticmethod
    def _lat_attr_setter(lat, attr_row):
        """
        sets new lattice parameters based attr_row
        and by thant updates the related PDFgenerator that is linked to 'lat'
        """
        for lat_attr, lat_val in attr_row.items():
            setattr(lat, lat_attr, lat_val)

    @staticmethod
    def _lin_combo(components, weights):
        """ Creates vectors of the linear combination of the components """
        combo = np.zeros(len(components.T))
        for comp, weight in zip(components, weights):
            combo += comp * weight
        return combo

    def _convolve_to_voigt(self, intensity):
        """
        convolve peaks to voigt

        Parameters
        ----------
        intensity: array
            xrd intensity

        Returns
        -------
        convolved intensity array
        """
        t = np.arange(-1000, 1000, 0.1)
        fwhmg = 1.5
        fwhml = 1.5
        v = sps.voigt_profile(t, fwhmg, fwhml)
        return np.convolve(intensity, v, 'same')

    def _get_mix_lat_attr_strains_df(self):
        mix_lat_attr_strains = dict()
        for lat_attr, i_max_strain in zip(self.lat_attrs, self.avg_max_strain):
            if i_max_strain != 0:
                mix_lat_attr_strains.update(
                    {lat_attr: np.linspace(0, i_max_strain, num=self.steps, endpoint=True).round(3)})
            else:
                mix_lat_attr_strains.update(
                    {lat_attr: np.zeros(self.steps).round(3)})
        df = pd.DataFrame(mix_lat_attr_strains)
        return df

    def _plot(self, datadict):
        myplot = Plot()
        for compname in self.comps_to_sink:
            for label, yvector in datadict[compname].items():
                myplot.plot(self.xvector, yvector, f'{compname}: {label}', num=compname)
        myplot.show()

    # =========================== CORE ====================

    def get_base_xrds(self):
        print("GENERATING XRD")
        cm = CifMap()
        for comp in self.compounds.items():
            print(f'generating: {comp}')
            compname = comp[0]  # name of the generator / phase
            print(f'calculating: "{compname}"')
            file = eval(f'self.{compname}')  # path to structure
            stru = loadCrystal(file)  # get structure
            gen = self.set_generator(compname, stru)
            self.set_gen_config(gen)  # set general config
            lat_df = self.get_lat_vals(gen)

            # for XRD
            dd = cm.load_cif(file, type='dd')
            dp = cm.load_cif(file, type='dp')
            cm.set_cell_dp2dd(dd, dp)
            config = dict(uiso=self.uiso)
            dd = cm.set_atoms(dd, dp, config=config)
            # nanop_envelop = self.generate_nano_particle_envelop(compname)
            new_uiso = self.uiso
            _uiso = list()
            for i in tqdm(range(self.steps)):
                _uiso.append(new_uiso)
                new_uiso += self.uiso_step
                dd.Atoms.uiso = new_uiso  # FIXME  not sure this is the right syntax
                # set new lattice dimensions
                attr_row = lat_df.iloc[i]
                self._lat_attr_setter(dd.Cell, attr_row)
                # simulate
                dd.Scatter.setup_scatter('xray')
                q, I = dd.Scatter.generate_powder(q_max=self.qmax, peak_width=0, background=0, powder_average=True)
                I = self._convolve_to_voigt(I)
                # archive
                self.xrddict[compname].update({str(attr_row.to_dict()): I})
        print("used Uiso:", _uiso, '\n')

        # setup x
        self.xvector = q

    def get_base_pdfs(self):
        print("GENERATING PDF")
        self.xvector = np.arange(self.rmin, self.rmax, self.rstep)
        for comp in self.compounds.items():
            print(f'generating: {comp}')
            compname = comp[0]  # name of the generator / phase
            print(f'calculating: "{compname}"')
            file = eval(f'self.{compname}')  # path to structure
            stru = loadCrystal(file)  # get structure
            gen = self.set_generator(compname, stru)
            self.set_gen_config(gen)  # set general config
            lat_df = self.get_lat_vals(gen)
            nanop_envelop = self.generate_nano_particle_envelop(compname)
            new_uiso = self.uiso
            _uiso = list()
            for i in tqdm(range(self.steps), leave=True, ascii=True):
                _uiso.append(new_uiso)
                self.set_uiso(gen, new_uiso)
                lat = gen.phase.getLattice()  # set new lattice dimensions
                attr_row = lat_df.iloc[i]
                self._lat_attr_setter(lat, attr_row)
                pdf = gen(self.xvector)
                if nanop_envelop:
                    pdf = pdf * nanop_envelop
                self.pdfdict[compname].update({str(attr_row.to_dict()): pdf})
                new_uiso += self.uiso_step
            print("used Uiso:", _uiso, '\n')

    def add_mixture_into_pdfdict(self):
        self.compounds.update(self.mix_compounds)
        self.pdfdict.update({self.mixname: dict()})
        similarity_counter = 1
        _used_mixtures = list()
        mix_lat_attr_strains_df = self._get_mix_lat_attr_strains_df()
        _weights = list()
        for i in range(self.steps):
            mix_lat_attr_strains_row = mix_lat_attr_strains_df.iloc[i]
            # strain_percentage = np.linspace(0, self.avg_max_strain, num=self.steps, endpoint=True) * 100
            components = list()
            for comp in self.compounds.items():
                compname = comp[0]
                if compname is not self.mixname:
                    arrays = self.pdfdict[compname]
                    array = list(arrays.values())[i]
                    components.append(array)
            components = np.asarray(components)
            combo_pdf = self._lin_combo(components, self.weights)
            mixture_name = mix_lat_attr_strains_row.to_dict()
            mixture_name.update({"w": np.array(self.weights).round(4).tolist()})
            if mixture_name in _used_mixtures:
                mixture_name = f'{mixture_name}_{similarity_counter}'
                similarity_counter += 1
            _used_mixtures.append(mixture_name)
            self.pdfdict[self.mixname].update({str(mixture_name): combo_pdf})
            # finish
            if self.SHIFT_WEIGHTS:
                self.weights = self.weights + self.weight_deltas
                _weights.append(list(self.weights))
        print("used WEIGHTS: ", _weights)

    def add_mixture_into_xrddict(self):
        self.compounds.update(self.mix_compounds)
        self.xrddict.update({self.mixname: dict()})
        similarity_counter = 1
        _used_mixtures = list()
        mix_lat_attr_strains_df = self._get_mix_lat_attr_strains_df()
        for i in range(self.steps):
            mix_lat_attr_strains_row = mix_lat_attr_strains_df.iloc[i]
            # strain_percentage = np.linspace(0, self.avg_max_strain, num=self.steps, endpoint=True) * 100
            components = list()
            for comp in self.compounds.items():
                compname = comp[0]
                if compname is not self.mixname:
                    arrays = self.xrddict[compname]
                    array = list(arrays.values())[i]
                    components.append(array)
            components = np.asarray(components)
            combo_xrd = self._lin_combo(components, self.weights)
            mixture_name = mix_lat_attr_strains_row.to_dict()
            mixture_name.update({"w": np.array(self.weights).round(4).tolist()})
            if mixture_name in _used_mixtures:
                mixture_name = f'{mixture_name}_{similarity_counter}'
                similarity_counter += 1
            _used_mixtures.append(mixture_name)
            self.xrddict[self.mixname].update({str(mixture_name): combo_xrd})
            # finish
            if self.SHIFT_WEIGHTS:
                print(self.weights)
                self.weights = self.weights + self.weight_deltas

    def to_csv(self, datadict, force_overwrite: bool = False, type: str = ''):
        dump = False
        for compname, compvals in datadict.items():
            if compname in self.comps_to_sink:
                comp_df = pd.DataFrame(compvals)
                comp_df.insert(0, 'xvector', self.xvector, True)
                filename = "_".join(map(str, [type, compname, self.lat_attrs, self.steps,
                                              self.compounds[compname]['max_strain']])) + ".csv"
                dirname = compname.split('{', 1)[0]
                dir = os.path.join(self.dumpdir, dirname)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
                dump_path = os.path.join(self.dumpdir, dir, filename)
                if os.path.isfile(dump_path):
                    while True:
                        if force_overwrite:
                            answer = 'y'
                        else:
                            answer = input(f'"{dump_path}" File exists.'
                                           f'\nOverwrite? [y/n]')
                        if answer.lower() == 'y':
                            dump = True
                            break
                        elif answer.lower() == 'n':
                            break
                        else:
                            continue
                else:
                    dump = True
                if dump:
                    comp_df.to_csv(dump_path)

    def to_txt(self, datadict, force_overwrite: bool = False, type: str = ''):
        dump = False
        for compname, compvals in datadict.items():
            if compname in self.comps_to_sink:
                comp_df = pd.DataFrame(compvals)
                comp_df.insert(0, 'xvector', self.xvector, True)
                filename = "_".join(
                    map(str, [type, compname, f'size[A]:{self.psize[compname]}', self.lat_attrs, self.steps,
                              self.compounds[compname]['max_strain']])) + ".npt"
                dirname = compname.split('{', 1)[0]
                dir = os.path.join(self.dumpdir, dirname)
                if not os.path.isdir(dir):
                    os.mkdir(dir)
                dump_path = os.path.join(self.dumpdir, dir, filename)
                if os.path.isfile(dump_path):
                    while True:
                        if force_overwrite:
                            answer = 'y'
                        else:
                            answer = input(f'"{dump_path}" File exists.'
                                           f'\nOverwrite? [y/n]')
                        if answer.lower() == 'y':
                            dump = True
                            break
                        elif answer.lower() == 'n':
                            break
                        else:
                            continue
                else:
                    dump = True
                if dump:
                    header = f"[{type}]\n"
                    header += "[compounds]\n"
                    header += yaml.dump(self.compounds, sort_keys=False, indent=2)
                    header += "\n[mixed compound]\n"
                    header += yaml.dump(self.mix_compounds, sort_keys=False, indent=2)
                    header += "\n[simulation configs]\n"
                    header += yaml.dump(self.sim_config, sort_keys=False, indent=2)
                    header += "\n[columns]\n"
                    header += yaml.dump(comp_df.columns.tolist(), sort_keys=False, indent=2)
                    np.savetxt(dump_path, comp_df.to_numpy().T, newline='\n\n',
                               fmt='%.4e', delimiter=', ', header=header)

    def main(self):

        # simulate PDFs or XRDs for each component
        if self.PDF:
            # simulate each component
            self.get_base_pdfs()
            # add linear combo
            self.add_mixture_into_pdfdict()
            self.combos_df = pd.DataFrame(self.pdfdict)
            # plot
            if self.PLOT:
                self._plot(self.pdfdict)
            # dump
            if self.DUMP:
                if self.DUMP_CSV:
                    self.to_csv(self.pdfdict, force_overwrite=self.FORCE_OVERWRITE, type='PDF')
                if self.DUMP_TXT:
                    self.to_txt(self.pdfdict, force_overwrite=self.FORCE_OVERWRITE, type='PDF')
        if self.XRD:
            # simulate each component
            self.get_base_xrds()
            # add linear combo
            self.add_mixture_into_xrddict()
            self.combos_df = pd.DataFrame(self.xrddict)
            # plot
            if self.PLOT:
                self._plot(self.xrddict)
            # dump
            if self.DUMP:
                if self.DUMP_CSV:
                    self.to_csv(self.xrddict, force_overwrite=self.FORCE_OVERWRITE, type='XRD')
                if self.DUMP_TXT:
                    self.to_txt(self.xrddict, force_overwrite=self.FORCE_OVERWRITE, type='XRD')


### IGNORE -- dev
def __demo_nanoparticle():
    from diffpy.srfit.pdf.characteristicfunctions import sphericalCF
    from matplotlib import pyplot as plt
    envelop = list()
    for r in Sim().xvector:
        envelop.append(float(sphericalCF(r, 10)))
    nanop_envelop = np.array(envelop)
    bulk = list(Sim().pdfdict['ZnSe'].values())[0]
    plt.plot(bulk, label='bulk')
    plt.plot(bulk * nanop_envelop, label='nanop => Eq: nanop_envelop * bulk')
    plt.plot(nanop_envelop * bulk[397], label='nanop_envelop')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    # gen = Sim().main()
    Sim().main()
