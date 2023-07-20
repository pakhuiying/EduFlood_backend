import mpmath
import cv2
import urllib.request
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import os
import math

def get_tiles(url):
    urllib.request.urlretrieve(url,"tile.png")
    tile = np.asarray(Image.open("tile.png"))

    return tile

def get_flood_tiles(x,y,z,e):
    """ 
    :params x,y,z: x and y tiles, z=zoom level
    :param e (f)
    """
    url_flood = f'https://www.floodmap.net/getFMTile.ashx?x={x}&y={y}&z={z}&e={e}'
    return get_tiles(url_flood)

def tile_num_to_lat_lon(xtile,ytile, zoom):
    """
    :param xtile (int): corresponds to google slippy map tile number
    :param ytile (int): corresponds to google slippy map tile number
    returns the coordinate of the upper left (northwest most) point of the tile,
    and the coordinate of the lower right (south east most) point of tile
    in WGS84 datum
    """
    def num2deg(xtile, ytile, zoom):
      n = 1 << zoom #2**zoom
      lon_deg = xtile / n * 360.0 - 180.0
      lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
      lat_deg = math.degrees(lat_rad)
      return lat_deg, lon_deg
    
    UL = num2deg(xtile, ytile, zoom)
    LR = num2deg(xtile+1, ytile+1, zoom)
    return UL,LR

def get_tile(lat_deg,lon_deg,zoom):
   ''' A function to get the relevant tile from lat,lon,zoom)'''
   
   lat_rad = math.radians(lat_deg)
   n =  1 << zoom
   
   xtile = n * ((lon_deg + 180) / 360)
   ytile = n * (1 - (math.log(math.tan(lat_rad) + mpmath.sec(lat_rad)) / math.pi)) / 2
   return int(xtile),int(ytile),int(zoom)

def lambert_distance(lat1, lon1, lat2, lon2):
    """ returns distance in meters between 2 points"""
    f = 1/298.257
    la1uLambert = lat1*math.pi/180
    la2uLambert = lat2*math.pi/180
    if (abs(lat1)<90):
        la1uLambert = math.atan((1 - f)*math.tan(la1uLambert))
    
    if (abs(lat2)<90):
        la2uLambert = math.atan((1 - f)*math.tan(la2uLambert))
    
    la1u = lat1*math.pi/180
    lo1u = lon1*math.pi/180
    la2u = lat2*math.pi/180
    lo2u = lon2*math.pi/180

    # Lambert Method
    deltaLat = la2uLambert - la1uLambert 
    deltaLon = lo2u - lo1u
    a = math.sin(deltaLat/2)*math.sin(deltaLat/2) + math.cos(la1uLambert)*math.cos(la2uLambert) * math.sin(deltaLon/2) * math.sin(deltaLon/2)
    c = 2*math.atan2(a**(1/2), (1 - a)**(1/2))
    P = (la1uLambert + la2uLambert)/2
    Q = (la2uLambert - la1uLambert)/2
    X = (c - math.sin(c))*math.sin(P)*math.sin(P)*math.cos(Q)*math.cos(Q)/math.cos(c/2)/math.cos(c/2) 
    Y = (c + math.sin(c))*math.sin(Q)*math.sin(Q)*math.cos(P)*math.cos(P)/math.sin(c/2)/math.sin(c/2)
    result_meters = 6378.1*(c - f*(X + Y)/2)*1000
    return result_meters

class AmISafe:
    def __init__(self, lat_deg,lon_deg,zoom):
        self.lat_deg = lat_deg
        self.lon_deg = lon_deg
        self.zoom = zoom

    def get_res(self,nrow,ncol,x,y,z):
        """ returns resolution per pixels"""
        UL,LR = tile_num_to_lat_lon(x,y,z)
        xres_deg = abs((LR[1] - UL[1])/ncol)
        yres_deg = abs((UL[0] - LR[0])/nrow)
        # print(yres_deg,xres_deg)
        y_meters = lambert_distance(UL[0], UL[1], LR[0], UL[1])
        x_meters = lambert_distance(UL[0], UL[1], UL[0], LR[1])
        print(y_meters, x_meters)
        yres_meters = y_meters/nrow
        xres_meters = x_meters/ncol
        return (UL,LR), (yres_deg,xres_deg), (yres_meters,xres_meters)

    def get_current_coord(self,UL,xres_deg, yres_deg):
        top = (UL[0] - self.lat_deg)/yres_deg
        left = (self.lon_deg - UL[1])/xres_deg
        return int(top), int(left)
    
    def get_radius(self, distance_meters, res):
        """ get radius in terms of pixels"""
        return math.ceil(distance_meters/res)
    
    def get_tiles(self, x,y,z,e1 = 1, e2 = 2):
        tile_1 = get_flood_tiles(x,y,z,e=e1)
        tile_2 = get_flood_tiles(x,y,z,e=e2)
        tiles = [tile_1[:,:,-1]/255,tile_2[:,:,-1]/255]
        return [t.astype(np.uint8) for t in tiles]
    
    def get_circles(self, nrow,ncol, center_coord, r1_pixel, r2_pixel):
        x, y = center_coord
        circles = []
        for r in [r1_pixel, r2_pixel]:
            im = np.zeros((nrow,ncol))
            image = cv2.circle(im, (x,y), radius=r, color=1, thickness=-1)
            circles.append(image)

        return circles
    
    def check_intersection(self,e1 = 1, e2 = 2,r1_meters = 50, r2_meters = 100, plot = True):
        x,y,z = get_tile(self.lat_deg, self.lon_deg, self.zoom)
        tiles = self.get_tiles(x,y,z,e1, e2)
        nrow, ncol = tiles[0].shape[0], tiles[0].shape[1]
        (UL,LR), (yres_deg,xres_deg), (yres_meters,xres_meters) = self.get_res(nrow,ncol,x,y,z)
        top, left = self.get_current_coord(UL,xres_deg, yres_deg)
    
        r1_pixel = self.get_radius(r1_meters, yres_meters)
        r2_pixel = self.get_radius(r2_meters, yres_meters)
        circles = self.get_circles(nrow,ncol, (left,top), r1_pixel, r2_pixel)
        
        itsxn_dict = {'elevation1': dict(), 'elevation2': dict()}
        for ki, i in zip(['elevation1','elevation2'], range(len(tiles))):
            for kj, j in zip(['radius1','radius2'],range(len(circles))):
                itsxn = np.sum(tiles[i]*circles[j]) # intersection
                itsxn_dict[ki][kj] = str(itsxn > 1)

        if plot is True:
            fig, axes = plt.subplots(2,len(tiles),figsize=(15,15))
            for ki, i in zip(['elevation1','elevation2'], range(len(tiles))):
                for kj, j in zip(['radius1','radius2'],range(len(circles))):
                    im = axes[j,i].imshow(tiles[i]+circles[j],vmin=0,vmax=2) # intersection
                    fig.colorbar(im, ax = axes[j,i])
                    axes[j,i].set_title(f'{ki}_{kj}')

            plt.show()

        return itsxn_dict
    
if __name__ == "__main__":
    ais = AmISafe(lat_deg=-31.970513,lon_deg=115.778649,zoom=15)
    itsxn_dict = ais.check_intersection(e1 = 1, e2 = 2,r1_meters = 50, r2_meters = 100, plot = True)
    print(f'itersection: {itsxn_dict}')