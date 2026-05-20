from myqueue.workflow import run


def workflow():
    run(script='solutions1.py', tmax='1h', cores=1)
