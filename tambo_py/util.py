import cmath 
import math  
import numpy as np
import numba
from scipy.interpolate import griddata
from numba import (  
    complex64,
    complex128,
    float32,
    float64,
    int32,
    int64,
    uint32,
    uint64,
    guvectorize,
    jit,
)

def make_topograhy_figure(...):
    import matplotlib.pyplot as plt
    #fig = plt.figure(figsize=(3,3))

    ax = plt.axes(projection='3d')
    x,y = np.meshgrid(longs,lats);
    g = griddata((Long, Lat), Elev, (x, y), method='cubic');
    g_del = np.delete(g,np.s_[1:11],1);
    g_gel = np.delete(g_del,np.s_[990:1000],1);
    y_del = np.delete(y,np.s_[990:1000],1);
    x_del = np.delete(x,np.s_[990:1000],1);
    #print(g)
    ax.contour3D(x_del[10:990],y_del[10:990],g_del[10:990], 400, cmap='terrain')
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_zlabel("elevation (m)")
    ax.view_init(70, 290)

    plt.gcf().set_size_inches(12, 12)
    #plt.savefig("ContourColcaFromAbove", bbox_inches='tight')
        