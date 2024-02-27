# =========================== SWITCHBOARD ===========================
class UserConfig:
    def __init__(self):
        # requested simulation
        self.XRD = False  # if True - simulate XRD
        self.PDF = True  # if True - simulate PDF

        # visualize results as plots?
        self.PLOT = True  # if true - visualize as plots

        # dump plots?
        self.DUMP = True  # if True - dump to dumpster (self.experiment)
        self.DUMP_CSV = True  # ovewrite as csv
        self.DUMP_TXT = True  # overwrite as txt (xy format)
        self.FORCE_OVERWRITE = True  # overwrite existing files

        # ----- set current working directory
        self.parent_path = '../example_data'  # RECOMMENDED:  set outside the `diffpysim` directory

        # ----- locate structure files
        self.cif_bank_dirname = 'cif_bank'  # this is the location of the cif structure files

        # ----- set dumpster
        self.experiment = 'ran_stretchnmf'  # your files will be saved here under `parent_path`

        # number of simulations
        self.steps = 20  # number of steps between zero strain to max strain

        # manipulate weights
        self.NORM_WEIGHTS = True  # normalize initial weights
        self.SHIFT_WEIGHTS = True  # include changing weights

        # ----- simulation config
        # XRD + PDF simulation parameters
        self.qdamp = 0.03  # term that increase a decay envelop to the PDF (relevant for high r; typically 0.03)
        self.qbroad = 0.000  # PDF peak broadening from increased intensity noise at high Q
        self.qmin = 0.1  # the smallest Q value; typically 0.1
        self.qmax = 30  # the largest Q value; at synchrotrons typically 20-30, but can be essentially any floating num.

        # additional parameters for PDF to simulations
        self.rmin = 0  # min value to of PDF
        self.rmax = 120  # max value to PDF
        self.rstep = 0.01  # # step interval

        # ----- since reused: the PDFGenerator's lattice parameters
        self.lat_attrs = ('a', 'b', 'c')  # set lattice attributes. You can change to single parameter, such as `a`

        # physical material properties  - ADP
        self.uiso = 0.005  # lowest (coldest) ADP
        self.uiso_max = 0.01  # highest (hotest). Used only to calculate uiso_step and as metadata
        self.compounds = {
            'ZnSe': {       # set compound name as the key
                'cifname': "ZnSe_216_F-43m_c.cif",  # the specfic filename.cif of the compound in the cif_bank dir
                'weight': 0.5,  # the initial phase fraction (range: 0.0-1.0)
                'max_weight_shift': 0,  # the shift in phase fraction (positive or negative, range: 0.0-1.0)
                'max_strain': (0.1, 0.1, 0.1),  # the shift in strain (expansion) matching self.lat_attrs
                'lat_attrs': self.lat_attrs,  # the PDFGenerator's lattice parameters
                'size': 0.0  # particle size in \AA ; for bulk, set size = 0
            },

            'BaTiO3_c': {
                'cifname': "BaTiO3_221_Pm-3m_c.cif",
                'weight': 0.5,
                'max_weight_shift': 0,
                'max_strain': (0.2, 0.2, 0.2),
                'lat_attrs': self.lat_attrs,
                'size': 0.0
            }
        }
