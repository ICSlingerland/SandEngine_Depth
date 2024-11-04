import numpy as np
import pyproj

def convert_rd_to_wgs84(points_rd):
    # Define the RD New (Amersfoort) CRS
    proj_rd_new = pyproj.Proj(init="epsg:28992")

    # Define the WGS84 CRS (used for lat/lon)
    proj_wgs84 = pyproj.Proj(init="epsg:4326")

    # Convert RD New to WGS84 (Lat, Lon) and retain height
    points_wgs84 = [
        (*pyproj.transform(proj_rd_new, proj_wgs84, point[0], point[1]), point[2])
        for point in points_rd
    ]
    
    return points_wgs84



