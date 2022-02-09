class Soln:

    """Soln class: reads and manipulates a Soln information taken from a Sinex solution
    Sinex (Solution/Technique Independant Exchange format) : Sinex format handling
    Description of the Sinex format is available at:
    http://www.iers.org/documents/ac/sinex/sinex_v210_proposal.pdf
    Soln allow to account for discontinuity in sites positions
    """
    
    def __init__(self,name=None):
        self.periods={}

    def append_period(self,soln,period):
        self.periods[soln]=period
        