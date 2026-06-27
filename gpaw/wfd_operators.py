import numpy as np

import gpaw.cgpaw as cgpaw
from gpaw.fd_operators import FDOperator


class WeightedFDOperator(FDOperator):
    def __init__(self, operators):
        """Compound operator A with nonlocal weights.

        A = Sum_i weights_i * operators_i

        Arguments:
        operators -- List of FDOperators.

        The weights have to be set using the set_weights method.
        A is build from operators which are then discarded.
        """
        self.gd = operators[0].gd
        self.shape = tuple(self.gd.n_c)
        self.coef_ps = []
        self.offset_ps = []
        for op in operators:
            assert self.gd == op.gd
            assert operators[0].dtype == op.dtype
            assert not hasattr(op, 'nweights')
            if 0 in op.offset_p:
                assert op.offset_p[0] == 0
                self.offset_ps.append(
                    np.ascontiguousarray(op.offset_p.copy())
                )
                self.coef_ps.append(
                    np.ascontiguousarray(op.coef_p.copy())
                )
            else:
                self.offset_ps.append(
                    np.ascontiguousarray(np.hstack(([0], op.offset_p))))
                self.coef_ps.append(
                    np.ascontiguousarray(np.hstack(([.0], op.coef_p))))
            assert self.offset_ps[-1][0] == 0
            assert len(self.coef_ps[-1]) == len(self.offset_ps[-1])
            assert self.offset_ps[-1].flags.c_contiguous
            assert self.coef_ps[-1].flags.c_contiguous
            assert self.offset_ps[-1].dtype == int
            assert self.coef_ps[-1].dtype == float

        # Ensure weights reside on the correct computational backend device (CPU vs GPU)
        self.xp = operators[0].xp
        self.nweights = len(operators)
        self.cfd = all([op.cfd for op in operators])
        self.mp = max([op.mp for op in operators])
        self.dtype = operators[0].dtype
        if self.gd.comm.size > 1:
            self.comm = self.gd.comm.get_c_object()
        else:
            self.comm = None
        self.npoints = max([op.npoints for op in operators])
        self.weights = None
        self.operator = None
        self.description = 'Weighted Finite Difference Operator\n      '
        self.description += '\n      '.join([op.description
                                             for op in operators])

        self.xp = np

    def set_weights(self, weights):
        """Set the operator weights.

        weights   -- List of numpy/cupy arrays sized gd.n_c
                     (weights are not copied).
        """
        assert len(weights) == self.nweights
        for weight in weights:
            assert weight.shape == self.shape
            
            # Check type matching dynamically based on computational backend
            if self.xp is np:
                assert weight.dtype == float
                assert weight.flags.c_contiguous
            else:
                from gpaw.gpu import cupy as cp
                assert weight.dtype == cp.float64 if hasattr(cp, 'float64') else float
        self.weights = weights
        
        # When compiling C level structures, pass correct representation
        if self.xp is np:
            self.operator = cgpaw.WOperator(
                self.nweights, self.weights,
                self.coef_ps, self.offset_ps, self.gd.n_c, self.mp,
                self.gd.neighbor_cd, self.dtype == float, self.comm, self.cfd
            )
        else:
            # Under GPU executions, CuPy manages memory allocations dynamically
            self.operator = cgpaw.WOperator(
                self.nweights, [w.get() for w in self.weights],
                self.coef_ps, self.offset_ps, self.gd.n_c, self.mp,
                self.gd.neighbor_cd, True, self.comm, self.cfd
            )
