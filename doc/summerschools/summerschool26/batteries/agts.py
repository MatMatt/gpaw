import os
import shutil
from pathlib import Path
from myqueue.workflow import run


def workflow():
    if os.getenv('AGTS_FILES'):
        dir = Path(os.getenv('AGTS_FILES'))
        for file in [Path('lifepo4_wo_li.traj'),
                     Path('NEB_init.traj')]:
            if not file.is_file():
                shutil.copyfile(dir / file, file)

    # batteries1:
    run(script='relax-graphite.py', tmax='1h')
    with run(script='solution3.py', tmax='1h', cores=8):
        run(script='solution4.py', tmax='1h', cores=8)
        run(script='solution5.py', tmax='1h', cores=8)

    # batteries2:
    d1 = run(script='fepo4.py', tmax='1h', cores=8)
    d2 = run(script='lifepo4.py', tmax='1h', cores=8)
    d3 = run(script='li_metal.py', tmax='1h', cores=8)
    run(script='eq_pot.py', deps=[d1, d2, d3])

    # batteries3:
    run(script='batteries3.py', tmax='1h', cores=8, deps=[s21, s22])
