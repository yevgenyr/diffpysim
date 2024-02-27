# this module responsible for mapping between different CIF readers.
# originated to map diffpy.structure.loadStructure(f) to pip install Dans-Diffraction.Crystal(f)

# xtl = dif.Crystal(f)
import Dans_Diffraction as ddif
import diffpy.structure
import os


class CifMap:
    """
    maps diffpy.structure onto Dan_Diffraction.Crystal
    -- name index --
        dd = Dan_diffraction
        dp = diffpy
    """

    def __init__(self):
        self.dirpath = ''
        self.dp_xtl = None
        self.dd_xtl = None

    def set_dirpath(self, dirpath: str):
        """set main dirpath"""
        assert os.path.exists(dirpath)
        self.dirpath = dirpath
        os.chdir(dirpath)

    @staticmethod
    def load_cif(filename: str, type: str = 'dd'):
        """
        Load structure

        Parameters
        ----------
        filename: the name of the cif file
        type: str
            type of loader. currently supports only dd and dp
            Default: 'dd'

        Returns
        -------

        """
        assert type in ['dd', 'dp']
        if type == 'dd':
            return ddif.Crystal(filename)
        elif type == 'dp':
            return diffpy.structure.loadStructure(filename)

    @staticmethod
    def set_cell_dp2dd(dd, dp):
        """
        set lattice parameters from dp to dd

        Parameters
        ----------
        dd: dd object
        dp: dp object

        Returns
        -------

        """
        dd.Cell.latt(dp.lattice.abcABG())

    def set_atoms(self, dd, dp, config:dict=None):
        """
        sets atom identity and properties from dd

        Parameters
        ----------
        dd: dd object
        dp: dp object
        config (optional): dict
            **conigs to dd.new_atoms

        Returns
        -------

        """
        _config = dict(u=dp.xyz.T[0],
                       v=dp.xyz.T[1],
                       w=dp.xyz.T[2],
                       type=dp.element,
                       label=dp.label,
                       occupancy=dp.occupancy,
                       uiso=dp.Uisoequiv,
                       mxmymz=None)
        if config:
            _config.update(config)
        dd.new_atoms(**_config)
        return dd

