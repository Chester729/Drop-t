class chromosomes:
    def __init__(self):
        self.human = [str(i) for i in range(1,23)] + ['X','Y']
        self.fly = ['2L','2R','2LHet','2RHet','3L','3R','3LHet','3RHet','4dm','U','Uextra','Xdm']
    def add_chr(self):
        return [f'chr{i}' for i in self]
    def remove_chr(self):
        return [i[3:] for i in self]