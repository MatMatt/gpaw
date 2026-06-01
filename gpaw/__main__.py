import gpaw
from gpaw.cli.main import main

# Set MPI backend to 'cgpaw' by default if executing
# `gpaw` or `python -m gpaw`
if gpaw.GPAW_MPI_BACKEND is None:
    gpaw.GPAW_MPI_BACKEND = "cgpaw"

if __name__ == '__main__':
    main()
