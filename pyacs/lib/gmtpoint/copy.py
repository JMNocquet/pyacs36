"""Copy method for GMT_Point."""


def copy(self):
    """Return a copy of the current GMT_Point."""
    NEW = self.__class__()
    NEW.code = self.code
    NEW.Ve = self.Ve
    NEW.Vn = self.Vn
    NEW.Vu = self.Vu
    NEW.SVe = self.SVe
    NEW.SVn = self.SVn
    NEW.SVen = self.SVen
    NEW.SVu = self.SVu
    NEW.lon = self.lon
    NEW.lat = self.lat
    NEW.he = self.he
    NEW.Cv_xyz = self.Cv_xyz
    NEW.Cv_enu = self.Cv_enu
    NEW.index = self.index
    return NEW
