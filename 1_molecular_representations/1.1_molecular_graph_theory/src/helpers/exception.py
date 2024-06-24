class NonexistentNodeError(Exception):
    def __init__(self, message: str = "Node does not exist in graph"):
        self.message = message
        super().__init__(self.message)


class NonexistentEdgeError(Exception):
    def __init__(self, message: str = "Edge does not exist in graph"):
        self.message = message
        super().__init__(self.message)


class CycleError(Exception):
    def __init__(self, message: str = "Cycle detected in graph"):
        self.message = message
        super().__init__(self.message)


class InvalidSmilesError(Exception):
    def __init__(self, message: str = "SMILES string is invalid"):
        self.message = message
        super().__init__(self.message)
