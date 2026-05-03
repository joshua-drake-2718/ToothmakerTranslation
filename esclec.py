from coreop2d import Coreop2d
import numpy as np

class Esclec():
    map: int
    read_failed: int
    ma_max: int = 5000
    positionsp = np.zeros((1000, 3, ma_max)) # attention if num_active_cells>1000 the system crashes when reading
    parap = np.zeros((32, ma_max))
    param_names = np.empty(32, dtype=object)
    knotsp = np.zeros((1000, ma_max), dtype=np.int32)
    ma: np.ndarray
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
        # Parameter file format (1-based for human readability):
        #   lines 1..30 map to parap[2..31, map].
        #   parap[0..1, map] are reserved for temps and num_active_cells,
        #   set elsewhere (set_params reads them but they default to 0).
        cls.read_failed = 0
        try:
            with open(filepath, 'r') as f:
                for i in range(2, 32):  # 30 slots: indices 2..31
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
    def set_params(cls, core: Coreop2d, imap: int):
        core.temps =            int(cls.parap[0, imap])  # global counter of simulation time steps
        core.num_active_cells = int(cls.parap[1, imap])  # current number of cells
        core.Egr = cls.parap[2, imap]
        core.Mgr = cls.parap[3, imap]
        core.Rep = cls.parap[4, imap]
        # parap[5, imap] reserved (corresponds to FORTRAN parap(6) — no field in the param file)
        core.Adh = cls.parap[6, imap]
        core.Act = cls.parap[7, imap]
        core.Inh = cls.parap[8, imap]
        # parap[9, imap] reserved
        core.Sec = cls.parap[10, imap]
        # parap[11, imap] reserved
        core.Da  = cls.parap[12, imap]
        core.Di  = cls.parap[13, imap]
        core.Ds  = cls.parap[14, imap]
        # parap[15, imap] reserved
        core.Int = cls.parap[16, imap]
        core.Set = cls.parap[17, imap]
        core.Boy = cls.parap[18, imap]
        core.Dff = cls.parap[19, imap]
        core.Bgr = cls.parap[20, imap]
        core.Abi = cls.parap[21, imap]
        core.Pbi = cls.parap[22, imap]
        core.Bbi = cls.parap[23, imap]
        core.Lbi = cls.parap[24, imap]
        core.Rad = int(cls.parap[25, imap])
        core.Deg = cls.parap[26, imap]
        core.Dgr = cls.parap[27, imap]
        core.Ntr = cls.parap[28, imap]
        core.Bwi = cls.parap[29, imap]
        core.ina = cls.parap[30, imap]
        core.umgr = cls.parap[31, imap]

    @classmethod
    def initialize_from_parameter_file(cls, core, filepath):
        cls.map = 0
        cls.read_param_file(filepath)
        if cls.read_failed == 0:
            cls.set_params(core, cls.map)

    @classmethod
    def guardaveinsoff_2(cls, core: Coreop2d, temp_neigh: np.ndarray, f=None):
        c = np.zeros(4)

        cls.ma = np.zeros(core.num_active_cells)
        cls.mat(core)

        def count_or_emit_faces(emit: bool):
            """Mirror of FORTRAN guardaveinsoff_2 face-emission loops (13.f90:1619-1656 / 1665-1721)."""
            nfaces = 0
            for i in range(core.num_active_cells):
                for j in range(core.nv_max):
                    bi = 0
                    nop = -1
                    ii = core.neigh[i, j]
                    if ii == -1 or ii >= core.num_active_cells:
                        continue
                    # ele: triangle search via i -> ii -> iii -> i
                    cycle_ale = False
                    for k in range(core.nv_max):
                        iii = core.neigh[ii, k]
                        if iii == -1 or iii >= core.num_active_cells or iii == i:
                            continue
                        cycle_ele = False
                        for kk in range(core.nv_max):
                            iiii = core.neigh[iii, kk]
                            if iiii == -1 or iiii >= core.num_active_cells:
                                continue
                            if iiii == i:    # triangle found
                                nfaces += 1
                                bi += 1
                                nop = iii
                                if emit and f is not None:
                                    f.write(f"3 {i} {ii} {iii}\n")
                                if bi == 1:
                                    cycle_ele = True
                                    break  # cycle ele
                                else:
                                    cycle_ale = True
                                    break  # cycle ale
                        if cycle_ale:
                            break
                        if cycle_ele:
                            continue
                    if cycle_ale:
                        continue
                    # quad search (only if not jumped to next j)
                    found_quad = False
                    for k in range(core.nv_max):
                        iii = core.neigh[ii, k]
                        if iii == -1 or iii >= core.num_active_cells or iii == i or iii == nop:
                            continue
                        if bi == 0:
                            continue
                        for kk in range(core.nv_max):
                            iiii = core.neigh[iii, kk]
                            if iiii == -1 or iiii >= core.num_active_cells or iiii == ii or iiii == nop:
                                continue
                            for kkk in range(core.nv_max):
                                jj = core.neigh[iiii, kkk]
                                if jj == i:    # quad found
                                    nfaces += 1
                                    if emit and f is not None:
                                        f.write(f"4 {i} {ii} {iii} {iiii}\n")
                                    found_quad = True
                                    break
                            if found_quad: break
                        if found_quad: break
            return nfaces

        # --- Pass 1: count faces ---
        nfaces = count_or_emit_faces(emit=False)

        # --- Write header and vertices ---
        if f is not None:
            f.write("COFF\n")
            f.write(f"{core.num_active_cells} {nfaces} 0\n")

            for i in range(core.num_active_cells):
                cls.get_rainbow_knot(core, cls.ma[i], cls.va_min, cls.va_max, c, i)
                x = core.positions[i, 0]
                y = core.positions[i, 1]
                z = core.positions[i, 2]
                f.write(f"{x} {y} {z} {c[0]} {c[1]} {c[2]} {c[3]}\n")

        # --- Pass 2: write faces ---
        count_or_emit_faces(emit=True)

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
            c[0], c[1], c[2], c[3] = 0.6, 0.6, 0.6, 0.8
        elif f < 0.2:
            c[0], c[1], c[2], c[3] = 1.0, f, 0.0, 0.5
        elif f < 1.0:
            c[0], c[1], c[2], c[3] = 1.0, f * 3, 0.0, 1.0
        else:
            c[0], c[1], c[2], c[3] = 1.0, 1.0, 0.0, 1.0

    @classmethod
    def get_rainbow_knot(cls, core, val, min_val, max_val, c, i):
        if max_val > min_val:
            f = (val - min_val) / (max_val - min_val)
        else:
            f = 0.5

        if core.knots[i] == 1:
            c[0], c[1], c[2], c[3] = 0, 1, 1, 0
        elif f < 0.07:
            c[0], c[1], c[2], c[3] = 0.6, 0.6, 0.6, 0.8
        elif f < 0.2:
            c[0], c[1], c[2], c[3] = 1.0, f, 0.0, 0.5
        elif f < 1.0:
            c[0], c[1], c[2], c[3] = 1.0, f * 3, 0.0, 1.0
        else:
            c[0], c[1], c[2], c[3] = 1.0, 1.0, 0.0, 1.0