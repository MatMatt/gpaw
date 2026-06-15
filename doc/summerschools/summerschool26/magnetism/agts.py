from myqueue.workflow import run


def workflow():
    with run(script='xyz.py'):
        run(script='CrI3_gs.py')
        run(script='CrI3_fm.py')
        run(script='get_Tc_mf.py')
        run(script='get_Tc_fit.py')
        run(script='CrI3_anisotropy.py')
        run(script='VI2_gs.py')
        run(script='VI2_afm.py')
        run(script='VI2_noncol.py')
        run(script='VI2_anisotropy.py')
