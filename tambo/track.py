from os import set_inheritable
from numba.core.decorators import njit
import numpy as np
from numba import jit

from geometry import Direction, Geometry, Point

class Track: 
    """
    Primary Particle Trajectory TAMBO.
    """  
    def __init__(self, point: Point, direction: Direction): 
        """
        Constructs particle trajectory.
        """    
        self.start_track = self.__get_start_track()
        self.point = self.start_track[0]
        self.direction = self.start_track[1]

    def __get_start_track(self): 
        """
        Builds 6 item array defining starting position and direction of track
        """    
        line_eq = [[self.point.x,self.point.y,self.point.z],
                   [self.direction.x_dir,self.direction.y_dir,self.direction.z_dir]]
        return np.around(np.array(line_eq, dtype = np.float32),3)
    
    def find_end_points(self, geometry: Geometry):
        """
        Returns the end points of the track for a given geometry.
        """   

        for i,boundaries in enumerate(geometry.geometry_box):

            start_point = self.start_track[0][i%3]
            start_direct = self.start_track[1][i%3]
        
            flag = False
        
            if start_direct == 0: 
                continue       
        
            t = (boundaries - start_point)/start_direct
        
            if t < 0 or t == 0:  
                continue

            pot_end = ([self.point.x+(t*self.point.x_dir),self.point.y+(t*self.point.y_dir),
                                self.point.z+(t*self.point.z_dir)])
            pot_end = np.around(np.array(pot_end,dtype = np.float32),3)
            
            for j,points in enumerate(pot_end): 

                if points < geometry.geometry_box[2*j] or points > geometry.geometry_box[2*j+1]:
                    flag = True 
                    break 
                
            if flag == False: 
                data = np.around(np.array(pot_end, dtype=np.float32),3)
                end_pts = np.around(np.subtract(data,self.start_track[0]),3)
                return self.start_track[0],end_pts

if __name__ == "__main__":
    print("building geometry")
    geo = Geometry("../resources/ColcaValleyData.txt")
    point = geo.coordinate_points[0]
    dir = Direction(90,90)
    print("building track")
    track = Track(point,dir)
    print(track.find_end_points(geo))