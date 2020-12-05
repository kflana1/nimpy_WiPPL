import numpy as np
from scipy.interpolate import griddata

'''
this module deals with the nimpy grid and geometry conversions
it depends on:
    numpy (as np)
    scipy.interpolate.griddata

in nimpy three geometries are allowed:
    cylindrical: (R, Phi, Z) cylindrical coordinates
    spherical: (R, Phi, Theta) spherical coordinates (theta is defined as |/ angle)
    cartesian: (X, Y, Z)

for nimrod, R&Z or R&Theta or X&Z are the poloidal components and 
Phi or Y is the toroidal component
'''
def convert_field_to_circ(theta, fR, fZ):
    '''
    Converts a field from rectangular to spherical components
    input: theta (np.ndarray dim=1) theta array to be used, needs to be same length as one dim of f
           fR (np.ndarray dim=2) R component of field
           fZ (np.ndarray dim=2) Z component of field
    output: Ar (np.ndarray dim=2) new spherical R component of field
            ATheta (np.ndarray dim=2) new theta component of field
    '''
    if fR.shape != theta.shape:
        diff = len(fR.shape) - len(theta.shape)
        if diff == 0:
            pass
        elif diff == 1:
            theta = theta[:,:,np.newaxis]*np.ones_like(fR)
        elif diff == 2:
            theta = theta[:,:,np.newaxis,np.newaxis]*np.ones_like(fR)
    AR = fR * np.sin(theta) + fZ * np.cos(theta)
    ATheta = fR * np.cos(theta) - fZ * np.sin(theta)
    return AR, ATheta 

def convert_field_to_rect(theta, fR, fTheta):
    '''
    Convers a field from spherical to rectangular components
    input: theta (np.ndarray dim=1) theta array to be used, needs to be same length as one dim of f
           fR (np.ndarray dim=2) spherical R component of field
           fTheta (np.ndarray dim=2) Theta component of field
    output: AR (np.ndarray dim=2) new cylindrical R component of field
            AZ (np.ndarray dim=2) new Z component of field
    '''
    if fR.shape != theta.shape:
        diff = len(fR.shape) - len(theta.shape)
        if diff == 1:
            theta = theta[:, :, np.newaxis]*np.ones_like(fR)
        elif diff == 2:
            theta = theta[:, :, np.newaxis, np.newaxis]*np.ones_like(fR)
    AR = fR * np.sin(theta) + fTheta * np.cos(theta)
    AZ = fR * np.cos(theta) - fTheta * np.sin(theta)
    return AR, AZ

def convert_grid_to_spherical(R, Z):
    ''' convert a R,Z grid to R,Theta

    Args:
        R (np.ndarray): cylindrical R or cartesian X
        Z (np.ndarray): cylindrical or cartesian Z

    Returns:
        R, Theta (tuple of np.ndarray's): 
            shperical_R=sqrt(R^2+Z^2)
            Theta = arctan(R/Z)
    '''
    return np.sqrt(R**2+Z**2) , np.arctan2(R,Z)

def convert_grid_to_rect(R, Theta):
    ''' convert a R, Theta grid to R, Z

    Args:
        R (np.ndarray): spherical R
        Theta (np.ndarray): spherical Theta

    Returns:
        R, Z (tuple of np.ndarray's):
            cylindrical_R=R*sin(Theta)
            Z=R*cos(Theta)
    '''
    return R*np.sin(Theta), R*np.cos(Theta)

def convert_grid(old_grid, old_geom, new_geom):
    ''' converts a grid to a new geometry

    converts a grid to a new geometry based on geometry names

    Args:
        old_grid (tuple of np.ndarray's): (R,Z), (R,Theta) or (X,Z) grid
        old_geom (str): 'cylindrical', 'spherical', or 'cartesian'
        new_geom (str): geom to change into

    Returns:
        (R, Z) or (R, Theta) or (X, Z) based on new_geom
    '''    
    if old_geom == 'cylindrical':
        if new_geom == 'cylindrical':
            print('already in cylindrical geometry')
            new_grid = old_grid
        elif new_geom == 'spherical':
            new_grid = convert_grid_to_spherical(old_grid[0], old_grid[1])
        elif new_geom == 'cartesian':
            new_grid = old_grid
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    elif old_geom == 'spherical':
        if new_geom == 'cylindrical' or new_geom == 'cartesian':
            new_grid = convert_grid_to_rect(old_grid[0], old_grid[1])
        elif new_geom == 'spherical':
            print('already in spherical geometry')
            new_grid = old_grid
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    elif old_geom == 'cartesian':
        if new_geom == 'cylindrical':
            new_grid = old_grid
        elif new_geom == 'spherical':
            new_grid = convert_grid_to_spherical(old_grid[0], old_grid[1])
        elif new_geom == 'cartesian':
            print('already in cartesian geometry')
            new_grid = old_grid
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    else:
        raise ValueError('{0} is not a valid old_geom. Try cylindrical, spherical or cartesian'.format(old_geom))
    return new_grid

def convert_field(f1_old, f3_old, old_grid, old_geom, new_geom):
    ''' converts a vector field to a new geometry

    converts the poloidal components of a vector field to a new
    geometry based on the geometry names. toroidal component is 
    the same for all geometries

    Args:
        f1_old (np.ndarray): R or X field component
        f3_old (np.ndarray): Z or Theta field component
        old_grid (tuple): (R, Z), (R, Theta) or (X,Z) grid
        old_geom (str): 'cylindrical', 'spherical', or 'cartesian'
        new_geom (str): geom to change into

    Returns:
        (f1, f3) tuple:
            f1 (np.ndarray): R or X field component
            f3 (np.ndarray): Z or Theta field component
    '''        
    if old_geom == 'cylindrical':
        if new_geom == 'cylindrical':
            print('already in cylindrical geometry')
            f1 = f1_old
            f3 = f3_old
        elif new_geom == 'spherical':
            _, Theta = convert_grid_to_spherical(old_grid[0], old_grid[1])
            f1, f3 = convert_field_to_circ(Theta, f1_old, f3_old) 
        elif new_geom == 'cartesian':
            f1 = f1_old
            f3 = f3_old
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    elif old_geom == 'spherical':
        if new_geom == 'cylindrical' or new_geom == 'cartesian':
            f1, f3 = convert_field_to_rect(old_grid[1], f1_old, f3_old)
        elif new_geom == 'spherical':
            print('already in spherical geometry')
            f1 = f1_old
            f3 = f3_old
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    elif old_geom == 'cartesian':
        if new_geom == 'cylindrical':
            f1 = f1_old
            f3 = f3_old
        elif new_geom == 'spherical':
            _, Theta = convert_grid_to_spherical(old_grid[0], old_grid[1])
            f1, f3 = convert_field_to_circ(Theta, f1_old, f3_old)
        elif new_geom == 'cartesian':
            print('already in cartesian geometry')
            f1 = f1_old
            f3 = f3_old
        else:
            raise ValueError('{0} is not a valid new_geom. Try cylindrical, spherical or cartesian'.format(new_geom))
    else:
        raise ValueError('{0} is not a valid old_geom. Try cylindrical, spherical or cartesian'.format(old_geom))
    return f1, f3

def interp(f, R, Z, new_grid, method='linear'):
    '''
    Interpolates data from a field onto a linear array, such as a probe chord.
    Input: f (np.ndarray dim=2) input field
           R (np.ndarray dim=2) cylindrical R grid values of f
           Z (np.ndarray dim=2) Z grid values of f
           newgrid (tuple of np.ndarray) new grid values (R,Z) to interpolate to
                    output field will be the same shape
           method (optional, default='linear') scipy.interpolate.griddata interpolation method
    '''
    points = list(zip(R.flatten(), Z.flatten()))
    return griddata(points, f.flatten(), new_grid, method=method)

def convert_vec_fields_to_circ(in_dict):
    '''
    Takes a dictionary of fields and transforms vector fields to appropriate grid coordinates.
    This will return the values on the grid output by nimplot, which can be unstructured. If you are having trouble plotting,
    try interpolating.
    dependencies: numpy (import numpy as np); interp (defined below)
    
    :param in_dict: (dict) input dictionary, must have a grid R, Z 
    :return: (dict) output dictionary with fields cast in spherical geometry and 'circ_R' and 'circ_Theta' coordinates
    '''

    out_dict = {}
    r,theta = convert_grid_to_spherical(in_dict['R'],in_dict['Z'])
    out_dict['circ_R'] = r
    out_dict['circ_Theta'] = theta
    vec_keys = [x for x in list(in_dict.keys()) if '_R' in x or '_Z' in x or '_Phi' in x]
    scalar_keys = [x for x in list(in_dict.keys()) if x not in vec_keys]
    vec_keys = list(set([x.replace('_R','').replace('_Z','').replace('_Phi','') for x in vec_keys]))
    for k in vec_keys:
        out_dict[k+'_R'], out_dict[k+'_Theta'] = convert_field_to_circ(theta, in_dict[k+'_R'], in_dict[k+'_Z'])
        out_dict[k+'_Phi'] = in_dict[k+'_Phi']
    for k in scalar_keys:
        out_dict[k] = in_dict[k]

    return out_dict


def make_regular_grid(in_dict, n=None, method='linear'):
    '''
    Takes a dictionary produced by read_vec (and stripped), interpolates it onto either a circular or rectangular
    grid, and transforms vector fields to appropriate grid coordinates.
    dependencies: numpy (import numpy as np); interp (defined below)
    
    :param in_dict: (dict) input dictionary, must have 'R' and 'Z' keys at the very least 
    :param n: (list; default=None) number of grid points to interpolate to. If n is None
                    the grid will be a square rounded to have approximately the same number of grid points as the original. 
                    If len(n)=1 will make a square grid with that number
                    of points. If len(n)=2 will make a rectangular grid with nR=n[0] and nZ=n[1].
    :return: (dict) updated dictionary cast onto regular grid.
    
    '''
    if n is None:
        n = 2*[int(np.sqrt(in_dict['R'].shape))]
    elif type(n) is not list:
        n = 2*[int(n)]
    if len(n) < 2:
        n = 2*[n[0]]
    out_dict = {}
    r = np.linspace(in_dict['R'].min(), in_dict['R'].max(), n[0])
    z = np.linspace(in_dict['Z'].min(), in_dict['Z'].max(), n[1])
    R, Z = np.meshgrid(r, z)
    points = list(zip(in_dict['R'], in_dict['Z']))
    keys = [kl for kl in list(in_dict.keys()) if kl not in ['R', 'Z']]
    for k in keys:
        if np.any(in_dict[k]):
            out_dict[k] = griddata(points, in_dict[k], (R, Z), method=method)
        else:
            out_dict[k] = np.zeros_like(R)
    out_dict['R'] = R
    out_dict['Z'] = Z

    return out_dict
