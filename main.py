import sys
import argparse
import os
from indexing import *
from coreop2d import Coreop2d
from esclec import Esclec


def main():
    parser = argparse.ArgumentParser(
        description='ToothMaker tooth morphogenesis simulator (Python translation)'
    )
    parser.add_argument('input_file', help='Parameter file path')
    parser.add_argument('output_folder', help='Folder for output files')
    parser.add_argument('output_name', help='Base name for output files')
    parser.add_argument('iterations', type=int,
                        help='Number of iterations per save block')
    parser.add_argument('save_blocks', type=int,
                        help='Number of save blocks')
    args = parser.parse_args()

    cac = args.input_file
    caufolder = args.output_folder
    cau = args.output_name
    iteration_total = args.iterations
    sstep = args.save_blocks

    core = Coreop2d
    io = Esclec

    core.max_z_layers = 4
    core.initial_conditions()
    io.initialize_from_parameter_file(core, cac)

    core.allocate_initial_state()
    prev_num_active_cells = core.num_active_cells
    io.set_params(core, 1)
    core.num_active_cells = prev_num_active_cells

    core.initact()

    os.makedirs(caufolder, exist_ok=True)

    for iti in range(1, abs(sstep) + 1):
        iter_label = str(iti * iteration_total)

        nff = os.path.join(caufolder, iter_label + '_' + cau)
        nff = nff.replace(' ', '_')

        nfoff = nff + '_.off'
        nfes = nff + '_.txt'

        core.temps = 0
        nt = iteration_total
        core.iteration(nt)

        with open(nfoff, 'w') as f:
            io.guardaveinsoff_2(core, core.neigh, f)

        print(f'Block {iti}/{abs(sstep)} complete: {nfoff}')


if __name__ == '__main__':
    main()
