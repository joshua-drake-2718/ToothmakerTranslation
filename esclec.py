from coreop2d import Coreop2d
from indexing import *

class Esclec():
    map: int
    read_failed: int
    ma_max: int = 5000
    positionsp = float_array((1000, 3, ma_max)) # attention if num_active_cells>1000 the system crashes when reading
    parap = float_array((32, ma_max))
    param_names = str_array(32)
    knotsp = int_array((1000, ma_max))
    ma: float_array
    va_max: float
    va_min: float
    cac: str
    cad: str
    caq: str
    cat: str
    cas: str
    cau: str
    cass: str

    @classmethod
    def read_param_file(cls):
        pass

    @classmethod
    def set_params(cls, core: Coreop2d, imap: int): # imap: index of parameter set (always 1 in current usage)
        core.temps = cls.parap[1, imap] # global counter of the number of simulation time steps that have elapsed, measured in number of iterations that have already been executed
        core.num_active_cells = cls.parap[2, imap] # Current number of cells
        core.Egr = cls.parap[3, imap]
        core.Mgr = cls.parap[4, imap]
        core.Rep = cls.parap[5, imap]
        # Nothing
        core.Adh = cls.parap[7, imap]
        core.Act = cls.parap[8, imap]
        core.Inh = cls.parap[9, imap]
        # Nothing
        core.Sec = cls.parap[11, imap]
        # Nothing
        core.Da = cls.parap[13, imap]
        core.Di = cls.parap[14, imap]
        core.Ds = cls.parap[15, imap]
        # Nothing
        core.Int = cls.parap[17, imap]
        core.Set = cls.parap[18, imap]
        core.Boy = cls.parap[19, imap] # Boy (Eq. 13)
        core.Dff = cls.parap[20, imap] # Dff (differentiation rate). (Eq. 6) (formerly tadif)
        core.Bgr = cls.parap[21, imap] # Bgr (border growth factor) (softeq. "Growth Biases", subroutine updatecellposition). (formerly fac)
        core.Abi = cls.parap[22, imap]
        core.Pbi = cls.parap[23, imap]
        core.Bbi = cls.parap[24, imap]
        core.Lbi = cls.parap[25, imap]
        core.Rad = cls.parap[26, imap]
        core.Deg = cls.parap[27, imap]
        core.Dgr = cls.parap[28, imap]
        core.Ntr = cls.parap[29, imap]
        core.Bwi = cls.parap[30, imap]
        core.ina = cls.parap[31, imap]
        core.umgr = cls.parap[32, imap]

    @classmethod
    def initialize_from_parameter_file(cls):
        pass

    @classmethod
    def guardaveinsoff_2(cls, core: Coreop2d, temp_neigh: int_array):
        c = float_array(4)
        mic = float_array(4)
        # I hace to go from neighbouring to faces

        face = int_array((core.num_active_cells*20, 5))
        nfaces: int
        bi: int
        nop: int

        cls.ma = float_array(core.num_active_cells)
        cls.mat()

        nfaces = 0

        # code inside ele loop
        # returns True if cycling ele
        def ele():
            nonlocal nfaces, bi, nop
            iii = core.neigh[ii, k]
            if iii == 0 or iii > core.num_active_cells or iii == i:
                return True # cycle ele
            for kk in range(core.nv_max):
                iiii = core.neigh[iii, kk]
                if iiii == 0 or iiii > core.num_active_cells:
                    return True # cycle ele
                if iiii == i: # triangle trobat
                    nfaces += 1
                    bi += 1
                    nop = iii
                    if bi == 1:
                        return True # cycle ele
                    return False # cycle ale

        for i in range(core.num_active_cells):
            for j in range(core.nv_max):
                bi = 0
                ii = core.neigh[i, j]
                if ii == 0 or ii > core.num_active_cells:
                    continue
                for k in range(core.nv_max):
                    if ele():
                        continue
                    else:
                        break # cycle ale
                else:
                    for k in range(core.nv_max):
                        iii = core.neigh[ii, k]
                        if iii == 0 or iii > core.num_active_cells or iii == i or iii == nop:
                            continue
                        if bi == 0:
                            continue
                        # to by the squares
                        for kk in range(core.nv_max):
                            iiii = core.neigh[iii, kk]
                            if iiii == 0 or iiii > core.num_active_cells or iiii == ii or iiii == nop:
                                continue
                            for kkk in range(core.nv_max):
                                jj = core.neigh[iiii, kkk]
                                if jj == i: # triangle found
                                    nfaces += 1
                                    break
                            else: continue
                            break
                        else: continue
                        break
        
        # TODO: write statements
        for i in range(core.num_active_cells):
            cls.get_rainbow_knot(cls.ma(i), cls.va_min, cls.va_max, c, i)

        nfaces = 0

        for i in range(core.num_active_cells):
            for j in range(core.nv_max):
                bi = 0
                ii = core.neigh[i, k]
                if ii == 0 or ii > core.num_active_cells:
                    continue
                for k in range(core.nv_max):
                    iii = core.neigh[ii, k]
                    if iii == 0 or iii > core.num_active_cells or iii == i:
                        continue
                    for kk in range(core.nv_max):
                        iiii = core.neigh[iii, kk]
                        if iiii == 0 or iiii > core.num_active_cells:
                            continue
                        if iiii == i: # triangle found
                            nfaces += 1
                            cls.get_rainbow(cls.ma[i], cls.va_min, cls.va_max, c)
                    


    @classmethod
    def mat(cls):
        pass

    @classmethod
    def get_rainbow(cls, val: float, min_val: float, max_val: float, c: float_array):
        pass

    @classmethod
    def get_rainbow_knot(cls, val: float, min_val: float, max_val: float, c: float_array, i: int):
        pass