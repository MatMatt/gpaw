from gpaw.dft import DFT, Parameters
from gpaw.new.pw.hybrids import ibz2bz
from gpaw.new.logger import Logger


def fixed(dft: DFT,
          bp,
          txt='-'):
    old_params = dft.params.todict()
    old_params.pop('h', None)
    kwargs = {**old_params,
              'gpts': dft.density.nt_sR.desc.size}
    params = Parameters(**kwargs)
    log = Logger(txt, dft.comm)
    builder = params.dft_component_builder(dft.atoms, log=log)
    basis_set = builder.create_basis_set()
    kbcomm1 = dft.ibzwfs.kpt_band_comm
    kbcomm2 = builder.communicators['D']
    potential = dft.potential.redist(
        builder.grid,
        builder.electrostatic_potential_desc,
        builder.atomdist,
        kbcomm1, kbcomm2)
    density = dft.density.redist(builder.grid,
                                 builder.interpolation_desc,
                                 builder.atomdist,
                                 kbcomm1, kbcomm2)
    ibzwfs = builder.create_ibz_wave_functions(basis_set, potential)
    # ibzwfs.fermi_levels = dft.ibzwfs.fermi_levels
    scf_loop = builder.create_scf_loop()
    scf_loop.update_density_and_potential = False
    scf_loop.fix_fermi_level = True  # not update_fermi_level
    for name in ['energy', 'density', 'forces']:
        scf_loop.convergence.pop(name, None)
    grid = dft.density.nt_sR.desc
    if 0:
        mypsits, _ = ibz2bz(
            dft.ibzwfs,
            dft.setups,
            dft.relpos_ac,
            grid=grid,
            plan=grid.fft_plans(),
            log=log)
    scf_loop.hamiltonian.update_wave_functions(dft.ibzwfs)
    dft = DFT.from_components(
        dft.atoms, ibzwfs, density, potential,
        builder.setups,
        scf_loop,
        builder.create_potential_calculator(),
        log,
        params=params,
        energies=dft.energies)

    return dft
