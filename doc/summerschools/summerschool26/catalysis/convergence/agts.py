from myqueue.workflow import run


def workflow():
    with run(script='eads.py', tmax='1h', cores=8)
    with run(script='convergence.py', tmax='1h', cores=8):
        run(script='plot.py')
