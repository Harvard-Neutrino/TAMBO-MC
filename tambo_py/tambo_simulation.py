from typing_extensions import ParamSpec
from numba.core.types.scalars import Float
import numpy as np
from scipy.integrate import quad

from tambo.geometry import Geometry, Direction, Point
from tambo.particle import Particle
from tambo.track import Track

class TAMBO(object):
    """
    Master class of TAMBO.
    """  
    def __init__(self, geometry: Geometry): 
        self.geometry = geometry

    def inject_particle(self, particle: Particle):
        end_points = particle.trajectory.find_end_points(self.geometry)
        column_density = self.__get_column_density(end_points)
        sampled_column_density = np.random.uniform()*column_density
        interaction_position = Point(sampled_column_density,...)
        particle.trajectory.point = interaction_position

    def __get_density(self, t:Float, array:np.ndarray):
    
        #array = beginning x,y,z and endpoint x,y,z 
    
        z = array[0][2] + (t * array[1][2])
        y = array[0][1] + (t * array[1][1])
        x = array[0][0] + (t * array[1][0])
    
        dzdt = array[1][2]
        dydt = array[1][1]
        dxdt = array[1][0]
    
        if z <= self.geometry.geometry_spline(x,y): 
            return self.geometry.density_of_rock * np.sqrt(1+ dzdt**2 + dydt**2 + dxdt**2)
        if z > self.geometry.geometry_spline(x,y): 
            return self.geometry.density_of_air * np.sqrt(1+ dzdt**2 + dydt**2 + dxdt**2)

    def __get_column_density(self,end_points:np.ndarray):
        return quad(lambda t: self.__get_density(t, end_points), 0, 1)

    def __get_interaction_point(self,column_density:Float)->Point:
        """
        Return interaction point.
        """     
        return Point(0.0,0.0,0.0)

    def __get_interaction_products(self,particle:Particle)->list(Particle):
        return

    def __inject_to_CORSIKA(self,particle_list:list(Particle)):
        return

if __name__ == "__main__":
    print("building geometry")
    geo = Geometry("../resources/ColcaValleyData.txt")
    point = geo.Coordinate_points[0]
    dir = Direction(0.1,0.1)
    print("building track")
    track = Track(point,dir)
    print("building TAMBO")
    tambito = TAMBO(geo)
    particle = Particle(0,10.0,track)
    print(tambito.inject_particle(particle))