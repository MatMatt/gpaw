from myqueue.workflow import run


def workflow():
    with run(script='check_convergence.py', tmax='1h', cores=8):
        run(script='plot.py')
