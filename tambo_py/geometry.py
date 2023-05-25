from os import set_inheritable
from numba.core.decorators import njit
from numba.core.types.iterators import EnumerateType
from numba.core.types.scalars import Float
import numpy as np
from numba import jit
from scipy.interpolate import SmoothBivariateSpline

#adding in jit for optimization 

class Point(object):
    """
    Coordinate point.
    """  
    def __init__(self, longitude: Float, latitude: Float, elevation: Float, lat_min: Float, long_min: Float):
        self.longitude = longitude
        self.latitude = latitude
        self.lat_min = lat_min 
        self.long_min = long_min 
        self.x, self.y = self.__set_distance_coordinate(latitude,longitude)
        self.z = elevation

#    @jit
    def __set_distance_coordinate(self, latitude: Float, longitude: Float):
        """
        Implements conversion from coordinate location to a coordinate in meters.
        """    
        lat_mid = (latitude + self.lat_min)/2.0
        m_per_deg_lat = (111132.954 - (559.822 * np.cos( 2.0 * lat_mid )) + (1.175 * np.cos( 4.0 * lat_mid)) + (0.0023 * np.cos( 6.0 * lat_mid)))
        m_per_deg_lon = (111412.82 * np.cos(lat_mid)) - (93.5*np.cos(lat_mid*3)) + (0.118*np.cos(5*lat_mid))
        delta_lat = latitude - self.lat_min 
        delta_long = longitude - self.long_min 

        x = delta_long * (m_per_deg_lon * 180/np.pi)
        y = delta_lat * (m_per_deg_lat * 180/np.pi)

        return x,y

class Direction(object):
    """
    Direction of motion.
    """    
    def __init__(self, phi: Float, theta: Float):
        self.phi = phi
        self.theta = theta

        self.x_dir = np.cos(theta)*np.sin(phi)
        self.y_dir = np.sin(theta)*np.sin(phi)
        self.z_dir = np.cos(phi)
        
class Geometry(object): 
    """
    Implements the basic geometry of TAMBO.
    """    
    def __init__(self,
                text, #text file of geometry data in lat/long/elev format 
                ): 

        datafile = np.genfromtxt(text, delimiter = "\t", skip_header = 1)

        self.lat = np.deg2rad(datafile[:,1])
        self.long = np.deg2rad(datafile[:,2])
        self.elev = datafile[:,3]
        self.__number_geo_points = len(self.lat)

        self.lat_max = np.max(self.lat)
        self.long_max = np.max(self.long)
        self.elev_max = np.max(self.elev)
        self.lat_min = np.min(self.lat)
        self.long_min = np.min(self.long)
        self.elev_min = np.min(self.elev)

        # A coord. system in meters with the origin at latmin, latmax
        self.coordinate_points = ([Point(self.long[i],self.lat[i],self.elev[i],self.lat_min, self.long_min) 
        for i in range(self.__number_geo_points)])

        self.geometry_spline = self.__construct_spline()
        self.geometry_box = self.__compute_dim_array()

        self.density_rock =  5520 #"kg/m^3"
        self.density_air = 1225 #"kg/m^3" 
        
#    @jit
    def __coords_to_meters(self, longitude: Float, latitude: Float): 
        """
        Implements conversion from coordinate location in degrees to a coordinate in meters.
        """    
        lat_mid = (latitude + self.lat_min)/2.0

        m_per_deg_lat = (111132.954 - (559.822 * np.cos( 2.0 * lat_mid )) + (1.175 * np.cos( 4.0 * lat_mid)) 
        + (0.0023 * np.cos( 6.0 * lat_mid)))
        m_per_deg_lon = (111412.82 * np.cos(lat_mid)) - (93.5*np.cos(lat_mid*3)) + (0.118*np.cos(5*lat_mid))
    
        delta_lat = latitude - self.lat_min 
        delta_long = longitude - self.long_min 
    
        x = delta_long * (m_per_deg_lon * 180/np.pi)
        y = delta_lat * (m_per_deg_lat * 180/np.pi)
        
        return np.around(np.array([x,y],dtype =np.float32),3)
    
    def __compute_dim_array(self): 
        """
        Returns a 6 item array defining the boundaries of geometry
        """    
        max_meters = self.__coords_to_meters(self.long_max,self.lat_max)
        array = ([0.0,max_meters[0],0.0,max_meters[1],self.elev_min,self.elev_max])
        return np.around(np.array(array,dtype =np.float32),3)
        
    def __construct_spline(self): 
        """
        Constructs a spline of geometry 
        """    
        x = [self.__coords_to_meters(g,i) for i,g in zip(self.lat,self.long)]
        meters_lat = [x[i][1] for i in range(len(x))]
        meters_long = [x[i][0] for i in range(len(x))]
        return SmoothBivariateSpline(meters_long,meters_lat,self.elev,kx=3,ky=3,s=100000)

if __name__ == "__main__":
    geo = Geometry("../resources/ColcaValleyData.txt")
    point = geo.coordinate_points[100]
    print(geo.geometry_spline(point.x,point.y),point.z)
    print(geo.geometry_box)