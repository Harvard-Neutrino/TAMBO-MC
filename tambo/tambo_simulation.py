#from typing_extensions import ParamSpec
from numba.core.types.scalars import Float
import numpy as np
from scipy.integrate import quad
import math 

from geometry import Geometry, Direction, Point
from particle import Particle
from track import Track

class TAMBO(object):
    """
    Master class of TAMBO.
    """  
    def __init__(self, geometry: Geometry): 
        self.geometry = geometry

    def inject_particle(self, particle: Particle):
        end_points = particle.trajectory.find_end_points(self.geometry)
        #return self.__get_column_density(end_points)
        column_density = self.__get_column_density(end_points)
        sampled_column_density = np.random.uniform()*column_density
        print(column_density, sampled_column_density)
        return self.__get_interaction_point(sampled_column_density,end_points)
        
        #interaction_position = Point(sampled_column_density,...)
        #particle.trajectory.point = interaction_position

    def __get_density(self, t:Float, array:np.ndarray):
    
        #array = beginning x,y,z and endpoint x,y,z 
    
        z = array[0][2] + (t * array[1][2])
        y = array[0][1] + (t * array[1][1])
        x = array[0][0] + (t * array[1][0])
    
        dzdt = array[1][2]
        dydt = array[1][1]
        dxdt = array[1][0]
    
        if z <= self.geometry.geometry_spline(x,y): 
            return self.geometry.density_rock * np.sqrt(1+ dzdt**2 + dydt**2 + dxdt**2)
        if z > self.geometry.geometry_spline(x,y): 
            return self.geometry.density_air * np.sqrt(1+ dzdt**2 + dydt**2 + dxdt**2)

    def __get_column_density(self,end_points:np.ndarray):
        integral, error = quad(lambda t: self.__get_density(t, end_points), 0, 1)
        return integral 

    def __get_interaction_point(self,column_density:Float,end_points)->Point:
        """
        Returns interaction point given a sampled column density
        """     
        length = np.linalg.norm(end_points[1])
        t = self.__binary_search(end_points,0,1,column_density,int(math.log(length,2)))

        z = end_points[0][2] + (t * end_points[1][2])
        y = end_points[0][1] + (t * end_points[1][1])
        x = end_points[0][0] + (t * end_points[1][0])

        #returns the coordinate in meters
        return x,y,z

    def __binary_search(self,end_points:np.ndarray,start,end,cd,step:int): 

        if end > start: 

            mid = (start+end)/2 

            if step == 0: 
                return mid 

            test_cd,error = quad(lambda t: self.__get_density(t, end_points), 0, mid)
            print(mid,test_cd,cd,step)
            if cd == test_cd: 
                return mid 
            if test_cd > cd: 
                return self.__binary_search(end_points,start,mid,cd,step-1)
            if test_cd < cd: 
                return self.__binary_search(end_points,mid,end,cd,step-1)

        else: 
            #something went wrong if this happens... 
            return -1 



    #def __get_interaction_products(self,particle:Particle)->list(Particle):
        #return

    #def __inject_to_CORSIKA(self,particle_list:list(Particle)):
        #return

if __name__ == "__main__":
    print("building geometry")
    geo = Geometry("../resources/ColcaValleyData.txt")
    point = geo.coordinate_points[9000]
    dir = Direction(90,45)
    print("building track")
    track = Track(point,dir)
    print("building TAMBO")
    tambito = TAMBO(geo)
    particle = Particle(0,10.0,track)
    print(tambito.inject_particle(particle))