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
    def read_param_file(cls, filepath):
        cls.read_failed = 0
        try:
            with open(filepath, 'r') as f:
                for i in range(3, 33):  # indices 3 to 32 inclusive
                    line = f.readline()
                    if not line:
                        cls.read_failed = 1
                        return
                    parts = line.split()
                    if len(parts) < 2:
                        cls.read_failed = 1
                        return
                    cls.parap[i, cls.map] = float(parts[0])
                    cls.param_names[i] = parts[1]
        except (IOError, ValueError):
            print("reading error for", filepath)
            cls.read_failed = 1

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
    def initialize_from_parameter_file(cls, core, filepath):
        cls.map = 1
        cls.read_param_file(filepath)
        if cls.read_failed == 0:
            cls.set_params(core, cls.map)

    @classmethod
    def guardaveinsoff_2(cls, core: Coreop2d, temp_neigh: int_array, f=None):
        c = float_array(4)
        mic = float_array(4)

        face = int_array((core.num_active_cells * 20, 5))
        nfaces = 0
        bi = 0
        nop = 0

        cls.ma = float_array(core.num_active_cells)
        cls.mat(core)

        # --- Pass 1: count faces ---
        def ele():
            nonlocal nfaces, bi, nop
            iii = core.neigh[ii, k]
            if iii == 0 or iii > core.num_active_cells or iii == i:
                return True
            for kk in range(core.nv_max):
                iiii = core.neigh[iii, kk]
                if iiii == 0 or iiii > core.num_active_cells:
                    return True
                if iiii == i:
                    nfaces += 1
                    bi += 1
                    nop = iii
                    if bi == 1:
                        return True
                    return False

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
                        break
                else:
                    for k in range(core.nv_max):
                        iii = core.neigh[ii, k]
                        if iii == 0 or iii > core.num_active_cells or iii == i or iii == nop:
                            continue
                        if bi == 0:
                            continue
                        for kk in range(core.nv_max):
                            iiii = core.neigh[iii, kk]
                            if iiii == 0 or iiii > core.num_active_cells or iiii == ii or iiii == nop:
                                continue
                            for kkk in range(core.nv_max):
                                jj = core.neigh[iiii, kkk]
                                if jj == i:
                                    nfaces += 1
                                    break
                            else: continue
                            break
                        else: continue
                        break

        # --- Write header and vertices ---
        if f is not None:
            f.write("COFF\n")
            f.write(f"{core.num_active_cells} {nfaces} 0\n")

            for i in range(core.num_active_cells):
                cls.get_rainbow_knot(core, cls.ma[i], cls.va_min, cls.va_max, c, i)
                x = core.positions[i, 1]
                y = core.positions[i, 2]
                z = core.positions[i, 3]
                f.write(f"{x} {y} {z} {c[1]} {c[2]} {c[3]} {c[4]}\n")

        # --- Pass 2: write faces ---
        nfaces = 0
        for i in range(core.num_active_cells):
            for j in range(core.nv_max):
                bi = 0
                ii = core.neigh[i, j]
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
                        if iiii == i:
                            nfaces += 1
                            bi += 1
                            nop = iii
                            if f is not None:
                                f.write(f"3 {i - 1} {ii - 1} {iii - 1}\n")
                            if bi == 1:
                                break
                            # quad face search
                            for k2 in range(core.nv_max):
                                iii2 = core.neigh[ii, k2]
                                if iii2 == 0 or iii2 > core.num_active_cells or iii2 == i or iii2 == nop:
                                    continue
                                if bi == 0:
                                    continue
                                for kk2 in range(core.nv_max):
                                    iiii2 = core.neigh[iii2, kk2]
                                    if iiii2 == 0 or iiii2 > core.num_active_cells or iiii2 == ii or iiii2 == nop:
                                        continue
                                    for kkk2 in range(core.nv_max):
                                        jj = core.neigh[iiii2, kkk2]
                                        if jj == i:
                                            nfaces += 1
                                            if f is not None:
                                                f.write(f"4 {i - 1} {ii - 1} {iii2 - 1} {iiii2 - 1}\n")
                                            break
                                    else: continue
                                    break
                                else: continue
                                break
                            break
                    else: continue
                    break

    @classmethod
    def mat(cls, core):
        cls.ma[:] = 0
        for i in range(core.num_active_cells):
            if core.knots[i] == 1:
                cls.ma[i] = 1.0
            else:
                if core.diff_state[i] > core.Int:
                    cls.ma[i] = 0.1
                if core.diff_state[i] > core.Set:
                    cls.ma[i] = 1.0
        cls.va_max = max(cls.ma[i] for i in range(core.num_active_cells))
        cls.va_min = min(cls.ma[i] for i in range(core.num_active_cells))

    @classmethod
    def get_rainbow(cls, val, min_val, max_val, c):
        if max_val > min_val:
            f = (val - min_val) / (max_val - min_val)
        else:
            f = 0.5

        if f < 0.07:
            c[1], c[2], c[3], c[4] = 0.6, 0.6, 0.6, 0.8
        elif f < 0.2:
            c[1], c[2], c[3], c[4] = 1.0, f, 0.0, 0.5
        elif f < 1.0:
            c[1], c[2], c[3], c[4] = 1.0, f * 3, 0.0, 1.0
        else:
            c[1], c[2], c[3], c[4] = 1.0, 1.0, 0.0, 1.0

    @classmethod
    def get_rainbow_knot(cls, core, val, min_val, max_val, c, i):
        if max_val > min_val:
            f = (val - min_val) / (max_val - min_val)
        else:
            f = 0.5

        if core.knots[i] == 1:
            c[1], c[2], c[3], c[4] = 0, 1, 1, 0
        elif f < 0.07:
            c[1], c[2], c[3], c[4] = 0.6, 0.6, 0.6, 0.8
        elif f < 0.2:
            c[1], c[2], c[3], c[4] = 1.0, f, 0.0, 0.5
        elif f < 1.0:
            c[1], c[2], c[3], c[4] = 1.0, f * 3, 0.0, 1.0
        else:
            c[1], c[2], c[3], c[4] = 1.0, 1.0, 0.0, 1.0