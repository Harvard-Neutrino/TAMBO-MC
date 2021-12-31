from os import set_inheritable
from numba.core.decorators import njit
from numba.core.types.iterators import EnumerateType
from numba.core.types.scalars import Float
import numpy as np
from numba import jit
from scipy.interpolate import SmoothBivariateSpline

class Point(object):
    """
    Coordinate point.
    """  
    def __init__(self,longitude: Float, latitude: Float, elevation: Float, latmin: Float, longmin: Float):
        self.Longitude = longitude
        self.Latitude = latitude
        self.X, self.Y = self.__set_distance_coordinate(latitude,longitude,latmin,longmin)
        self.Z = elevation

#    @jit
    def __set_distance_coordinate(self,latitude: Float,longitude: Float, latmin: Float, longmin: Float):
        """
        Implements conversion from coordinate location to a coordinate in meters.
        """    
        latMid = (latitude + latmin)/2.0
        m_per_deg_lon = (111132.954 - (559.822 * np.cos( 2.0 * latMid )) + (1.175 * np.cos( 4.0 * latMid)) + (0.0023 * np.cos( 6.0 * latMid)))
        m_per_deg_lat = (111412.82 * np.cos(latMid)) - (93.5*np.cos(latMid*3)) + (0.118*np.cos(5*latMid))
        delta_lat = latitude - latmin 
        delta_long = longitude - longmin 

        x = delta_long * (m_per_deg_lon * 180./np.pi)
        y = delta_lat * (m_per_deg_lat * 180./np.pi)

        return x,y

class Direction(object):
    """
    Direction of motion.
    """    
    def __init__(self,phi: Float,theta: Float):
        self.Phi = phi
        self.Theta = theta

        self.XDir = np.cos(theta)*np.sin(phi)
        self.YDir = np.sin(theta)*np.sin(phi)
        self.ZDir = np.cos(phi)
        
class Geometry(object): 
    """
    Implements the basic geometry of TAMBO.
    """    
    def __init__(self,
                text, #text file of geometry data in lat/long/elev format 
                ): 

        datafile = np.genfromtxt(text, delimiter = "\t", skip_header = 1)

        self.Lat = np.deg2rad(datafile[:,1])
        self.Long = np.deg2rad(datafile[:,2])
        self.Elev = datafile[:,3]
        self.__NumberGeoPoints = len(self.Lat)

        self.LatMax = np.max(self.Lat)
        self.LongMax = np.max(self.Long)
        self.ElevMax = np.max(self.Elev)
        self.LatMin = np.min(self.Lat)
        self.LongMin = np.min(self.Long)
        self.ElevMin = np.min(self.Elev)

        self.CoordinatePoints = [Point(self.Lat[i],self.Long[i],self.Elev[i],self.LatMin, self.LongMin) for i in range(self.__NumberGeoPoints)]

        self.GeometrySpline = self.__construct_spline()
        self.GeometryBox = self.__compute_dim_array()

        self.DensityRock =  5520 #"kg/m^3"
        self.DensityAir = 1225 #"kg/m^3" 
        
#    @jit
    def __coords_to_meters(self,longitude: Float,latitude: Float): 
        """
        Implements conversion from coordinate location in degrees to a coordinate in meters.
        """    
        latMid = (latitude + self.LatMin)/2.0

        m_per_deg_lat = (111132.954 - (559.822 * np.cos( 2.0 * latMid )) + (1.175 * np.cos( 4.0 * latMid)) 
        + (0.0023 * np.cos( 6.0 * latMid)))
        m_per_deg_lon = (111412.82 * np.cos(latMid)) - (93.5*np.cos(latMid*3)) + (0.118*np.cos(5*latMid))
    
        delta_lat = latitude - self.LatMin 
        delta_long = longitude - self.LongMin 
    
        x = delta_long * (m_per_deg_lon * 180/np.pi)
        y = delta_lat * (m_per_deg_lat * 180/np.pi)
        
        return np.around(np.array([x,y],dtype =np.float32),3)
    
    def __compute_dim_array(self): 
        """
        Returns a 6 item array defining the boundaries of geometry
        """    
        max_meters = self.__coords_to_meters(self.LongMax,self.LatMax)
        array = ([0.0,max_meters[0],0.0,max_meters[1],self.ElevMin,self.ElevMax])
        return np.around(np.array(array,dtype =np.float32),3)
        
    def __construct_spline(self): 
        """
        Constructs a spline of geometry 
        """    
        x = [self.__coords_to_meters(g,i) for i,g in zip(self.Lat,self.Long)]
        Meters_Lats = [x[i][1] for i in range(len(x))]
        Meters_Longs = [x[i][0] for i in range(len(x))]
        return SmoothBivariateSpline(Meters_Longs,Meters_Lats,self.Elev,kx=3,ky=3)

if __name__ == "__main__":
    geo = Geometry("../resources/ColcaValleyData.txt")
    point = geo.CoordinatePoints[1000]
    print(point.X)
    print(geo.GeometrySpline(point.X,point.Y),point.Z)
    print(geo.GeometryBox)