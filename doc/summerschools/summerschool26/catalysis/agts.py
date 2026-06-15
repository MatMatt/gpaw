# web-page: N2Ru_hollow.png, 2NadsRu.png, TS.xyz
from myqueue.workflow import run


def workflow():
    with run(script="check_convergence.py", tmax="1h", cores=8):
        run(script="convergence.py")

    with run(script="part1/n2_on_metal.py", tmax="2h"):
        with run(script="part3/neb.py", tmax="30h", cores=8):
            with run(script="part_extra/relax.py", tmax="12h", cores=24):
                run(script="part_extra/vibrations_initial.py", tmax="12h", cores=24)
                run(script="part_extra/vibrations_final.py", tmax="12h", cores=24)
            with run(
                script="part_extra/generate_transition_state_xyz.py", tmax="12h", cores=24
            ):
                run(script="part_extra/vibrations_transition.py", tmax="12h", cores=24)
