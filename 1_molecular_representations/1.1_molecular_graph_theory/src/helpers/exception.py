class CycleError(Exception):
    def __init__(self, message: str = "Cycle detected in graph"):
        self.message = message
        super().__init__(self.message)
