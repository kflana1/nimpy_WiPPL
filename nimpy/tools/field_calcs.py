import numpy as np
from scipy.integrate import cumtrapz
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

'''
this module contains all the vector calculus math used in nimpy on vector and scalar
fields. It depends on:
    numpy (definied as np)
    scipy.integrate.cumtrapz (as cumtrapz)
    warnings

Note: b/c 1/R blows up @ R=0, this would normally throw a lot of warnings, but 
they have been disabled.

----For all calculations d/dphi=0, i.e. we are only considering poloidal derivatives----

In nimpy, the R (or X) component is the axis=1 of numpy arrays and the Z component is the axis=0
Nimrod coordinate system has (R,Z) or (R, Theta) or (X, Z) as the poloidal components for cylindrical, spherical, and cartesian geometries, respectively. 
The Phi or Y component is the toroidal component for either cylindrical or spherical and cartesian geometries, respectively.
Break down:
    ----geom='cartesian'----
    f1 -> f_x
    f2 -> f_y "toroidal comp"
    f3 -> f_z
    ----geom='cylindrical'----
    f1 -> f_r
    f2 -> f_phi "toroidal comp"
    f3 -> f_z (same as f3 for cartesian)
    ----geom='spherical'----
    f1 -> f_r
    f2 -> f_phi "toroidal comp"
    f3 -> f_theta (polar angle)
'''

def df_dx1(f, x1):
    return np.gradient(f, axis=1)/np.gradient(x1, axis=1)

def df_dx3(f, x3):
    return np.gradient(f, axis=0)/np.gradient(x3, axis=0)

def d2f_dx12(f, x1):
    return np.gradient(df_dx1(f, x1), axis=1)/np.gradient(x1, axis=1)

def d2f_dx32(f, x3):
    return np.gradient(df_dx3(f, x3), axis=0)/np.gradient(x3, axis=0)

def div(f1, f3, x1, x3, geom='cylindrical'):
    '''
    Divergence of a nimpy field in poloidal plane
    input: f1 (np.ndarray) the R or X component of the field
           f3 (np.ndarray) the Z or Theta component of the field
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian
    '''
    if geom == 'cylindrical':
        return df_dx1(f1, x1) + f1/x1 + df_dx3(f3,x3)
    elif geom == 'spherical':
        return df_dx1(f1, x1) + (2./x1)*f1 + (df_dx3(f3, x3)/x1)
    elif geom == 'cartesian':
        return df_dx1(f1, x1) + df_dx3(f3, x3)
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def grad(f, x1, x3, geom='cylindrical'):
    '''
    Gradient of a nimpy field in poloidal plane
    input: f (np.ndarray) the field component
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian
    
    output based on geometry:
        cylindrical (dict) {'R': R component, 'Z': Z component}
        spherical (dict) {'R': R component, 'Theta': Theta component}
        cartesian (dict) {'X': X component, 'Z': Z component}
    '''
    if geom == 'cylindrical':
        return {'R': df_dx1(f, x1), 'Z': df_dx3(f, x3)}
    elif geom == 'spherical':
        return {'R': df_dx1(f, x1), 'Theta': df_dx3(f, x3)/x1}
    elif geom == 'cartesian':
        return {'X': df_dx1(f, x1), 'Z': df_dx3(f, x3)}
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def curl(f1, f2, f3, x1, x3, geom='cylindrical'):
    '''
    Curl of a nimpy field in poloidal plane
    input: f1 (np.ndarray) the R or X component of the field
           f2 (np.ndarray) the Phi or Y component of the field
           f3 (np.ndarray) the Z or Theta component of the field
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian

    output based on geometry:
        cylindrical (dict) {'R': R component, 'Phi': Phi component, 'Z': Z component}
        spherical (dict) {'R': R component, 'Phi': Phi component, 'Theta': Theta component}
        cartesian (dict) {'X': X component, 'Y': Y component, 'Z': Z component}
    '''
    if geom == 'cylindrical':
        ans = {}
        ans['R'] = -1.*df_dx3(f2, x3)
        ans['Phi'] = df_dx3(f1, x3)-df_dx1(f3, x1)
        ans['Z'] = df_dx1(f2, x1) + f2/x1
        return ans 
    elif geom == 'spherical':
        ans = {}
        ans['R'] = -1.*df_dx3(f2, x3)/x1
        ans['Phi'] = (df_dx3(f1,x3)/x1) - (f3/x1) - df_dx1(f3, x1) 
        ans['Theta'] = (1./x1)*f2 + df_dx1(f2, x1)
        return ans 
    elif geom == 'cartesian':
        ans = {}
        ans['X'] = -1.*df_dx3(f2, x3)
        ans['Y'] = df_dx3(f1, x3) - df_dx1(f3, x1)
        ans['Z'] = df_dx1(f2, x1)
        return ans 
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def laplacian(f, x1, x3, geom='cylindrical'):
    '''
    Laplacian of a nimpy field in poloidal plane
    input: f (np.ndarray) the field component
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian
    '''
    if geom == 'cylindrical':
        return d2f_dx12(f, x1) + df_dx1(f, x1)/x1 + d2f_dx32(f, x3)
    elif geom == 'spherical':
        return d2f_dx12(f, x1) + (2./x1)*df_dx1(f,x1) + (1./x1**2)*d2f_dx32(f, x3) 
    elif geom == 'cartesian':
        return d2f_dx12(f, x1) + d2f_dx32(f, x3)
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def vec_laplacian(f1, f2, f3, x1, x3, geom='cylindrical'):
    '''
    Vector Laplacian of a nimpy field in poloidal plane
    input: f1 (np.ndarray) the R or X component of the field
           f2 (np.ndarray) the Phi or Y component of the field
           f3 (np.ndarray) the Z or Theta component of the field
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian

    output based on geometry:
        cylindrical (dict) {'R': R component, 'Phi': Phi component, 'Z': Z component}
        spherical (dict) {'R': R component, 'Phi': Phi component, 'Theta': Theta component}
        cartesian (dict) {'X': X component, 'Y': Y component, 'Z': Z component}
    '''

    if geom == 'cylindrical':
        ans = {}
        ans['R'] = laplacian(f1, x1, x3) - f1/x1**2
        ans['Phi'] = laplacian(f2, x1, x3) - f2/x1**2
        ans['Z'] = laplacian(f3, x1, x3)
        return ans
    elif geom == 'spherical':
        ans = {}
        ans['R'] = laplacian(f1, x1, x3, geom=geom) - (2./x1**2)*(f1 + df_dx3(f3, x3))
        ans['Phi'] = laplacian(f2, x1, x3, geom=geom) - f2/x1**2 
        ans['Theta'] = laplacian(f3, x1, x3, geom=geom) - f3/x1**2 + (2./x1**2)*df_dx3(f1, x3)
        return ans 
    elif geom == 'cartesian':
        ans = {}
        ans['X'] = laplacian(f1, x1, x3, geom=geom)
        ans['Y'] = laplacian(f2, x1, x3, geom=geom)
        ans['Z'] = laplacian(f3, x1, x3, geom=geom)
        return ans 
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def f_dot_grad_f(f1, f2, f3, x1, x3, geom='cylindrical'):
    '''
    A_dot_grad_A of a nimpy field in poloidal plane
    input: f1 (np.ndarray) the R or X component of the field
           f2 (np.ndarray) the Phi or Y component of the field
           f3 (np.ndarray) the Z or Theta component of the field
           x1 (np.ndarray) the R grid
           x3 (np.ndarray) the Z or Theta grid
           geom (str): geometry to evaluate on, choices are: cylindrical, spherical or cartesian

    output based on geometry:
        cylindrical (dict) {'R': R component, 'Phi': Phi component, 'Z': Z component}
        spherical (dict) {'R': R component, 'Phi': Phi component, 'Theta': Theta component}
        cartesian (dict) {'X': X component, 'Y': Y component, 'Z': Z component}
    '''

    if geom == 'cylindrical':
        return {'R': f1*df_dx1(f1, x1) + f3*df_dx3(f1, x3) - f2**2/x1, 'Phi': f1*df_dx1(f2, x1) + f3*df_dx3(f2, x3) + (f1*f2)/x1, 'Z': f1*df_dx1(f3, x1) + f3*df_dx3(f3, x3)} 
    elif geom == 'spherical':
        ans = {}
        ans['R'] = f1*df_dx1(f1, x1) + (f3/x1)*df_dx3(f1, x3) - (f2**2 + f3**2)/x1
        ans['Phi'] = f1*df_dx1(f2, x1) + (f3/x1)*df_dx3(f2, x3) + (f2*f1)/x1
        ans['Theta'] = f1*df_dx1(f3, x1) + (f3/x1)*df_dx3(f3, x3) + (f3*f1)/x1
        return ans 
    elif geom == 'cartesian':
        ans = {}
        ans['X'] = f1*df_dx1(f1, x1) + f3*df_dx3(f1, x3)
        ans['Y'] = f1*df_dx1(f2, x1) + f3*df_dx3(f2, x3)
        ans['Z'] = f1*df_dx1(f3, x1) + f3*df_dx3(f3, x3)
        return ans 
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))

def calc_poloidal_stream_func(f3, x1, x3=None, geom='cylindrical'):
    '''
    Poloidal stream function calculator using scipy.integrate.cumtrapz
    input: f3 (np.ndarray) the Z or Theta component of the field
           x1 (np.ndarray) the R or X grid
           x3 (np.ndarray) (optional, default=None) if doing spherical calculation, needs the Theta grid
           geom (str): geometry to evalute on, choices are: cylindrical, spherical or cartesian
    '''
    if geom == 'cylindrical' or geom == 'cartesian':
        psi = cumtrapz(f3 * x1, x1[0, :], initial=0.0, axis=1)
    elif geom == 'spherical':
        if x3 is None:
            raise ValueError('You must provide a theta mesh for spherical stream calc')
        psi = -1.*np.sin(x3) * cumtrapz(f3 * x1, x1[0, :], initial=0.0, axis=1)
    else:
        raise ValueError('{0} is not a valid geom. Try cylindrical, spherical or cartesian'.format(geom))
    return psi
