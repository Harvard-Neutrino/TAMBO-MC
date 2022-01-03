from tambo.track import Track

class Particle(object): 
    """
    Particle Class representing a particle in TAMBO
    """
    def __init__(self, pdg_id: int, energy: float, trajectory:Track):
        self.pdg_id = pdg_id
        self.energy = energy
        self.trajectory = trajectory
        
if __name__ == "__main__":
    from tambo.geometry import Geometry, Direction, Point
    print("building track")
    dir = Direction(0.1,0.1)
    pon = Point(4,4,4,4,4)
    track = Track(pon,dir)
    print("building particle")
    particle = Particle(0,10.0,track)