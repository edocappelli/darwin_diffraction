class DiffractionException(Exception):
    def __init__(self, exception_text):
        super(DiffractionException, self).__init__(exception_text)


class ReflectionImpossibleException(DiffractionException):
    def __init__(self):
        super(ReflectionImpossibleException, self).__init__("Impossible geometry. Asymmetry angle larger than Bragg angle in Bragg geometry. No reflection possible.")