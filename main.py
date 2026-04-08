from indexing import *
from coreop2d import Coreop2d
from esclec import Esclec

sstep: int
nt: int
parapo = float_array(32)
cu: str
cd: str
ct: str
cq: str
acu: str
acd: str
act_str: str
acq: str
dcu: str
dcd: str
dct: str
dcq: str
nff: str
nfoff: str
nfes: str
accau_folder: str
cau_folder: str
dir_suffix: str
paras: int_array
sis: int
siss: int
idi: int

paras = int_array(buffer=[
    3, 4, 5, 7, 8,
    9, 11, 13, 14, 15,
    17, 18, 19, 20, 21, 23, 27, 28, 29
])