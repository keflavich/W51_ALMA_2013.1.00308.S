
def offset_to_point(xx, yy, path):
    """
    Determine the offset along the path to the nearest point on a path to the
    specified point
    """
    import shapely.geometry as geom

    if hasattr(path,'_coords'):
        line = geom.LineString(zip(path._coords.ra.deg, path._coords.dec.deg))
    else:
        raise ValueError("Path doesn't have coords - maybe need _xy?")
    point = geom.Point(xx, yy)
    return line.project(point)
