import gpaw
from gpaw.cli.main import main

# Set MPI backend to 'cgpaw' by default if executing
# `gpaw` or `python -m gpaw`
if gpaw.GPAW_MPI_BACKEND is None:
    import gpaw.cgpaw as cgpaw
    if cgpaw.have_mpi:
        gpaw.GPAW_MPI_BACKEND = "cgpaw"
    else:
        gpaw.GPAW_MPI_BACKEND = "serial"

if __name__ == '__main__':
    main()
