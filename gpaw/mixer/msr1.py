class MSR1Mixer(BaseMixer):
    def __init__(self,
                 nmaxold: int = 16,
                 beta: float = 0.08,
                 reg: float = 1e-4):
        self.nmaxold = nmaxold
        self.beta = beta
        self.reg = reg

    def _initialize_history_(self):
        ...

    def mix(self, Density: Density) -> float:
        ...
