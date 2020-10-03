import context

class Results:
    def __init__(self):
        self.source = None
        self.theta = None
        self.phi = None
        self.crystal = None
        self.SGlobal = None
        self.rx = None
        self.ry = None
        self.rz = None
        self.tx = None
        self.ty = None
        self.tz = None
        self.rTE = None # USEFUL PARAMETERS FOR ELLIPSOMETERY
        self.rTM = None
        self.tanPsi = None
        self.delta = None
        self.cosDelta = None
        self.R = None
        self.T = None
        self.RTot = None
        self.TTot = None
    def __str__(self):
        resultsString = f"Source ---\n \u03B8: {self.theta}, \u03C6: {self.phi}\n"
        resultsString += f"Reflection: rx: {self.rx:.3f}, ry: {self.ry:.3f} rz: {self.rz:.3f}\n"
        resultsString += f"Reflection: rTE: {self.rTE:.3f}, rTM: {self.rTM:.3f}\n"
        resultsString += f"RTot: {self.RTot:.3f}, TTot: {self.TTot:.3f}"
        return resultsString
