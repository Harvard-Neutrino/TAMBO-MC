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
        self.point = point 
        self.direction = direction
        self.start_track = self.__get_start_track()
        self.point_array = self.start_track[0]
        self.direction_array = self.start_track[1]

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

            #((i/2)+((-1)**(i)-1)/4) = 0,0,1,1,2,2
            start_point = self.start_track[0][int((i/2)+((-1)**(i)-1)/4)]
            start_direct = self.start_track[1][int((i/2)+((-1)**(i)-1)/4)]
        
            flag = False
        
            if start_direct == 0: 
                continue       
        
            t = (boundaries - start_point)/start_direct
            
            if t < 0 or t == 0:  
                continue

            pot_end = ([self.start_track[0][0]+(t*self.start_track[1][0]),
                        self.start_track[0][1]+(t*self.start_track[1][1]),
                        self.start_track[0][2]+(t*self.start_track[1][2])])

            pot_end = np.around(np.array(pot_end,dtype = np.float32),3)
            
            for j,potential_end in enumerate(pot_end): 

                if potential_end < geometry.geometry_box[2*j] or potential_end > geometry.geometry_box[2*j+1]:
                    flag = True 
                    break 
                
            if flag == False: 
                data = np.around(np.array(pot_end, dtype=np.float32),3)
                end_pts = np.around(np.subtract(data,self.start_track[0]),3)
                return self.start_track[0],end_pts

if __name__ == "__main__":
    print("building geometry")
    geo = Geometry("../resources/ColcaValleyData.txt")
    print(geo.geometry_box)
    point = geo.coordinate_points[9000]
    print(point.x,point.y)
    dir = Direction(90,45)
    print("building track")
    track = Track(point,dir)
    print(track.start_track)
    print("building track vector")
    print(track.find_end_points(geo))