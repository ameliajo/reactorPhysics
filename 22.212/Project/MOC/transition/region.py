import numpy as np

class Region():
    """The Region class describes the region using CSG

    Parameters
    ----------
    surfaces : list of Surface objects
        The surfaces that define the region
    orientations : list of -1,1
        On which side of each surface the region is located
    mat : int, opt
        id of material contained within surface
    name : str, opt
        name of the region

    Attributes
    ----------
    surfaces : list of Surface objects
        The surfaces that define the region
    orientations : list of -1,1
        On which side of each surface the region is located
    uid : int
        Unique id of the region
    mat : str, opt
        name of material contained within surface
    """
    def __init__(self, surfaces, uid, mat_class,mat = [],phi=0, orientations = [1, -1, -1, 1, 1],vol=0):
        self.surfaces = surfaces
        self.orientations = orientations
        self.uid = uid
        self.mat = mat
        self.phi = phi
        self.vol = vol
        self.mat_class = mat_class

        self.q = np.zeros([len(phi),])
        self.tracks_phi = np.zeros([len(phi),])
        self.q_phi = np.zeros([len(phi),])
        self.tot_track_length = 0

    def evaluate(cls, r):
        """Check if r is located within the region
        r : 2-tuple coordinates of point
        """
        for idx, surface in enumerate(cls.surfaces):
            sgn = cls.orientations[idx]
            if sgn*surface.evaluate(r) < 0:
                return False
        return True

def what_region(r, regions):
    """ Find out what region the point is in 
    Parameters
    ----------
    r : 2-tuple
        coordinates
    regions : list of Region objects
        The regions being used in this problem 

    Returns
    -------
    str 
        name of the region 
    """
    for region in regions:
        if region.evaluate(r):
            return region.uid
        
    
