import numpy as np
import matplotlib.pyplot as plt
from nimpy.tools import grid_geom, field_calcs, toroidal_modes, fileIO, plotting 

class NimData(object):
    ''' NimData class for manipulating data contained in nimpy h5 files

    The NimData object reads in scalar and vector fields from nimpy
    h5 files into a python object with many useful methods. This is
    the main class for any nimpy h5 file. It depends on the ScalarField
    and VectorField classes also contained in this module. 

    Attributes:
        h5name (str): the h5name used to create this object
        geom (str): the current geometry of the data in this object
            
            options are:
                'cylindrical' (default), 'spherical' and 'cartesian'
       
            based on geometry the following attributes are included:
            cylindrical:
                R (np.ndarray ndim=2): R values on poloidal grid
                Z (np.ndarray ndim=2): Z values on poloidal grid
                Phi (float, np.ndarray): Phi value(s) if not axisymmetric
            spherical:
                R (np.ndarray ndim=2): spherical R values on poloidal grid
                Theta (np.ndarray ndim=2): Theta values on poloidal grid
                Phi (float, np.ndarray): Phi value(s) if not axisymmetric
            cartesian (note: X=cylindrical R and Z=cylindrical Z):
                X (np.ndarray ndim=2): X values on poloidal grid
                Z (np.ndarray ndim=2): Y values on poloidal grid
                Y (float, np.ndarray): Y value(s) 

        comps (list of str): the names of the three vector components
        step (int): the time step of the data in this object
        time (float): the time in simulation units
        nmodes (int): number of toroidal modes 
        keff (np.ndarray): array of effective toroidal mode numbers
        inputs (dict): a dictionary of the nimset.out inputs used (only if -i
            option is used in nimh5 when making this h5 file)

        Fields:
           n (ScalarField): number density
           n0 (ScalarField): equilibrium number density
           P (ScalarField): total pressure
           P0 (ScalarField): equilibrium total pressure
           PE (ScalarField): electron pressure
           PE0 (ScalarField): equilibrium electron pressure
           tele (ScalarField): electron temperature
           tion (ScalarField): ion temperature
           diff (ScalarField): diff shape for elecd or isovisc params
           conc (ScalarField): material concentration (???)
           V (VectorField): velocity field
           V0 (VectorField): equilibrium velocity field
           B (VectorField): magentic field
           B0 (VectorField): equilibrium magnetic field
           J (VectorField): current density field
           J0 (VectorField): equilibrium current density field
    '''       
    def __init__(self, h5name, geom='cylindrical', phi=None):
        ''' NimData initialization function
        
        The init function sets the attributes for the NimData class based on 
        the h5 file name, the geometry and (if multiple modes are present) the 
        toroidal angle. 

        Args:
            h5name (str): the h5 filename containing nimpy data MUST INCLUDE EXTENSION
            geom (str): the geometry you wish to work with. Default is cylindrical
                because this corresponds to geom='tor' grid='rect' options in nimrod.in
                that WiPAL uses. Other options are 'spherical' and 'cartesian'. These
                geometries do not affect the grid, so plotting will be the same regardless.
                They only affect the attribute name space and the math methods in the
                ScalarField and VectorField classes
            phi (float, np.ndarray): the toroidal angle in radians that you wish to cast
                data onto. Only for cases with nmodes>1. The default angle is 0. If you 
                give a single float, the modes will be included on that plane. Else, if
                you provide an array, the modes will be cast onto those angles and the
                data will be 3D. The phi angle can be changed with the change_phi method.
        '''
        self.h5name = h5name
        big_dict = fileIO.load_h5_to_dict(h5name)
        ### if the input values are included in the h5 file, get them ###
        if 'input_vals' in list(big_dict.keys()):
            self.inputs = big_dict['input_vals']
        ### define the grid used for plotting ###
        self._grid = big_dict['grid']['R'], big_dict['grid']['Z']
        ### geometry specific assignments: comps, X1, X2, and X3 ###
        self._phi = phi
        if geom == 'cylindrical':
            self.comps = ['R', 'Phi', 'Z']
            self.R = self._grid[0]
            self.Z = self._grid[1]
            self.Phi = phi
        elif geom == 'spherical':
            self.comps = ['R', 'Phi', 'Theta']
            self.R, self.Theta = grid_geom.convert_grid_to_spherical(self._grid[0], self._grid[1])
            self.Phi = phi
        elif geom == 'cartesian':
            self.comps = ['X', 'Y', 'Z']
            self.X = self._grid[0]
            self.Z = self._grid[1]
            self.Y = phi
        else:
            raise ValueError('{0} is not a valid geometry. Try cylindrical, spherical or cartesian'.format(geom))
        self.geom = geom
        ### get the global dump values, but only make step, time, nmodes and keff public ###
        ### (others are nbl, nrbl, poly_deg, and nseams ###
        for k, f in list(big_dict['global_vals'].items()):
            if k in ['step', 'time', 'nmodes', 'keff']:
                setattr(self, k, f)
            else:
                setattr(self, '_'+k, f)
        ### this sets the toroidal angle to 0 if left on default and nmodes >1 ###
        if self.nmodes == 1:
            self.Phi = None
            self._phi = None
        elif phi is None:
            self._phi = 0.
            if geom == 'cylindrical' or 'spherical':
                self.Phi = 0.
            elif geom == 'cartesian':
                self.Y = 0.
        ### read in fields ###
        for field_type in ['eq_fields', 'fields']:
            for k, f in list(big_dict[field_type]['scalar'].items()):
                setattr(self, k, ScalarField(f, self._grid, self.geom, self.keff, self._phi))
            for vec_k, vec_f in list(big_dict[field_type]['vector'].items()):
                if self.geom is not 'cylindrical':
                    f1, f3 = grid_geom.convert_field(vec_f['R'], vec_f['Z'], self._grid, 'cylindrical', self.geom)
                else:
                    f1 = vec_f['R']
                    f3 = vec_f['Z']
                f2 = vec_f['Phi']
                setattr(self, vec_k, VectorField(f1, f2, f3, self._grid, self.geom, self.keff, self._phi))
                    
    
    def __str__(self):
        line1 = 'NimData object created from {0}'.format(self.h5name)
        if self.nmodes == 1:
            line2 = '\t'+'{0} geometry'.format(self.geom)
        elif self.geom in ['cylindrical', 'spherical']:
            phi_rad = self.Phi/np.pi
            if type(phi_rad) is np.ndarray:
                line2 = '{0} geometry with {1} toroidal modes @ phi=({2:.2f} - {3:.2f})*pi radians with {4} elements'.format(self.geom, self.nmodes, phi_rad[0], phi_rad[-1], phi_rad.size)
            else:
                line2 = '{0} geometry with {1} toroidal modes @ phi={2:.2f}*pi radians'.format(self.geom, self.nmodes, self.Phi)
        else:
            line2 = '\t'+'{0} geometry with {1} periodic modes @ Y={2:.2f}'.format(self.geom, self.nmodes, self.Y)
        line3 = '\n--------Global Values--------'+'\n'+'step: {0} \ntime: {1} \nnmodes: {2} \nkeff: {3}'.format(self.step, self.time, self.nmodes, self.keff)
        line4 = '\n--------Scalar Fields--------\n'
        for name, obj in list(self.__dict__.items()):
            if type(obj) is ScalarField:
                line4 += '{0}\n\t'.format(name)+obj.__str__()+'\n'
        line5 = '\n--------Vector Fields--------\n'
        for name, obj in list(self.__dict__.items()):
            if type(obj) is VectorField:
                line5 += '{0}\n\t'.format(name)+obj.__str__()+'\n'
        return line1+'\n'+line2+line3+line4+line5

    def __iter__(self):
        '''
        iter method will only iterate over the ScalarField and VectorField objects
        '''
        for attr, value in list(self.__dict__.items()):
            if type(value) in [ScalarField, VectorField]:
                yield attr, value

    def grid(self):
        ''' get the plotting grid

        Since nimpy saves on rect grids, this function will return a tuple of R, Z 
        that can be used for plotting. 

        Returns:
            R, Z (tuple of np.ndarray ndim=2)

        Example:
            R, Z = NimDataObj.grid()
            plt.contourf(R, Z, NimDataObj.tele.field())
        '''    
        return self._grid

    def tor_angle(self):
        ''' clearly prints the toroidal angle
        '''
        if self.nmodes == 1:
            print('axisymmetric')
        elif self.geom in ['cylindrical', 'spherical']:
            phi_rad = self.Phi / np.pi
            if type(phi_rad) in [np.ndarray, np.generic]:
                print(('phi= array from {0:.2f}*pi to {1:.2f}*pi radians with {2} elements'.format(phi_rad[0], phi_rad[-1], phi_rad.size)))
            else:
                print(('phi= {0:.2f}*pi radians'.format(phi_rad)))
        else:
            print(('not quite sure how this works in cartesian geometry...'+'\n'+'Y={0}'.format(self.Y)))

    def print_inputs(self):
        ''' clearly prints the input namelist if it was in h5 file

        nimh5 has an option (-i) to include input namelist used in a nimrod run.
        This method will clearly print this list (it is saved as a dictionary at
        attribute inputs)
        '''
        if hasattr(self, 'inputs'):
            fileIO.dic_print(self.inputs)
        else:
            print('No input dictionary included in this dump. Try running nimh5 again with -i option.')

    def add_ScalarField(self, name, f):
        ''' adds a new scalar field to this NimData object

        Allows you to add a new scalar field on the same grid as the
        other fields in this instance.

        Args:
            name (str): name of the new field (name of attribute to create)
            f (np.ndarray): new field on the same grid as current NimData instance
        '''
        setattr(self, name, ScalarField(f, self._grid, self.geom, self.keff, self._phi))

    def add_VectorField(self, name, f1, f2, f3):
        ''' adds a new vector field to this NimData object

        Allows you to add a new vector field on the same grid as the
        other fields in this instance.

        Args:
            name (str): name of the new field (name of attribute to create)
            f1 (np.ndarray): R, sphericalR, or X component of the field
            f2 (np.ndarray): Phi or Y component of the field
            f3 (np.ndarray): Z or Theta component of the field
        '''
        setattr(self, name, VectorField(f1, f2, f3, self._grid, self.geom, self.keff, self._phi))
    
    def recast_geom(self, new_geom):
        ''' recast the geometry of fields to a new one

        This method re-calls the initialization routine with a new geometry. It
        will update all of the fields' geometry in place. In order to generate a
        new instance with a different geometry, just re-call NimData.

        Args:
            new_geom (str): new geometry, must be 'cylindrical', 'spherical', or
                'cartesian'
        '''
        if new_geom == self.geom:
            print(('already in {0} geometry'.format(self.geom)))
        else:
            for name in list(self.__dict__.keys()):
                if name not in ['h5name', '_phi']:
                    delattr(self, name)
            self.__init__(self.h5name, geom=new_geom, phi=self._phi)

    def change_phi(self, new_phi):
        ''' change the toroidal angle for nmodes>1 data

        This method will change the poloidal plane to a new toroidal angle for 
        nmodes>1 data, or, if new_phi is a np.ndarray, it will cast the data into
        3D.

        Args:
            new_phi (float or np.ndarray): new toroidal angle in radians
        '''
        if self.nmodes == 1:
            print('disabled function for axisymmetric data')
            return
        for name, att in self.__iter__():
            att.__change_phi__(new_phi)
        self._phi = new_phi
        if self.geom in ['cylindrical', 'spherical']:
            self.Phi = new_phi
        else:
            self.Y = new_phi
    
class ScalarField(object):
    ''' ScalarField class for manipulating nimpy scalar fields

    This class allows for manipulation of scalar fields created by nimpy. It
    is an attribute of the NimData class for all instances of scalar fields.

    Attributes:
        field (np.ndarray): the field data. If nmodes>1 this is the sum of all
            of the mode data. 
        modes (np.ndarray): the mode data. Only exists if nmodes>1. This is the
            contribution from each mode indexed by axis=2
        min (float): the minimum value of the field
        max (float): the maximum value of the field
        avg (float): the mean value of the field
        nmodes (float): the number of toroidal modes for this field
        shape (tuple): the shape of the field attribute
    '''
    def __init__(self, field, grid, geom, keff, phi):
        ''' ScalarField initialization routine

        This method initializes an instance of ScalarField. It determines if there are
        more than one toroidal modes and then evaulates the complex field if nmodes>1.

        Args:
            field (np.ndarray): the field to initialize with
            grid (tuple of np.ndarray's): the plotting grid for the data
            geom (str): geometry to use for calculations. Can be 'cylindrical', 'spherical',
                or 'cartesian'
            keff (np.ndarray): the effective toroidal mode numbers (use [0] for axisymmetric)
            phi (float or np.ndarray): the toroidal angle(s) to cast the field to if nmodes>1
        '''
        self._grid = grid
        self._geom = geom
        if self._geom == 'spherical':
            self._coords = grid_geom.convert_grid_to_spherical(grid[0], grid[1])
        else:
            self._coords = self._grid
        self._keff = keff
        self._phi = phi
        self._org_field = field
        if keff.size > 1 and field.dtype not in [float, 'float16', 'float32', 'float64']:
            self.modes = toroidal_modes.eval_complex_field(field, keff, phi=phi)
            self.field = toroidal_modes.eval_complex_field(field, keff, phi=phi, sum=True)
        else:
            self.field = field
        self.min, self.max, self.avg, self.shape = self.__eval_field__(self.field)
    
    def __change_phi__(self, new_phi):
        self.__init__(self._org_field, self._grid, self._geom, self._keff, new_phi)

    def __eval_field__(self, field=None):
        ''' useful helper function for evaluating a field
        Args:
            field (np.ndarray): input field
        Returns:
            min, max, mean, shape (float, float, float, tuple)
        '''
        if field is None:
            field = self.field
        return np.nanmin(field), np.nanmax(field), np.nanmean(field), field.shape
        
    def __str__(self):
        line1 = 'ScalarField nimpy object, shape:{0}'.format(self.shape)
        if not np.any(self.field):
            line2 = 'empty array'
        else:
            line2 = 'min: {0}, max: {1}, avg: {2}'.format(self.min, self.max, self.avg)
        return line1+'\n\t'+line2

    def interpolate(self, new_grid, m=None, method='linear'):
        ''' interpolate field onto a new grid in poloidal plane

        This method can interpolate the field onto a new grid provided. If nmodes>1, it
        can interpolate either the sum of all the fields or an individual mode. Currently,
        this method is disabled for cases with phi as a np.ndarray. In these cases, change
        the phi to a single angle and then interpolate.

        Args:
            new_grid (tuple of nd.array): new grid points to interpolate onto. If in 
                cylindrical or cartesian geometry this is (R,Z) or (X,Z), else if in 
                spherical geometry this is (R, Theta). 
            m (int, default=None): If nmodes>1, this is the mode number to interpolate.
                If left to default, it will interpolate the sum of all modes. Doesn't
                do anything if nmodes=1.
            method (str, default='linear'): the interpolation method used by the call
                to scipy.interpolate.griddata
        
        Returns:
            interp_field (np.ndarray): the field interpolated onto the new grid
        '''
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
            ans = grid_geom.interp(self.field, self._coords[0], self._coords[1], new_grid, method=method)
        else:
            m = int(m)
            ans = grid_geom.interp(self.modes[:,:,m], self._coords[0], self._coords[1], new_grid, method=method)
        return ans
    
    def edge(self, m=None):
        ''' finds the edge values for a field

        This method steps through rows of Z looking for the last value before NaN. 

        Args:
            m (int, default=None): If nmodes>1, this is the mode number to interpolate.
                If left to default, it will interpolate the sum of all modes. Doesn't
                do anything if nmodes=1.
        
        Returns:
            r, edge_field (tuple of np.ndarrays): the r value of the edge and the field at the edge
                both will be the length of the Z array
        '''
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
           field = self.field
            
        else:
            m = int(m)
            field = self.modes[:,:,m]

        edge = np.zeros(self._grid[1].shape[0])
        r = np.zeros(self._grid[1].shape[0])
        for i in range(self._grid[1][:,0].size):
            row = self.field[i,:]
            ixs = np.where(np.isnan(row))[0]
            if ixs.size == 0:
                ix = -1
            else:
                ix = ixs[0] - 1
            r[i] = self._grid[0][i,ix]
            edge[i] = row[ix]
        
        return r,edge 

    def grad(self, m=None):
        ''' gradient of field in poloidal plane

        This method takes the gradient of the field in the correct geometry. Currently,
        this method is disabled for cases with phi as a np.ndarray.

        Args:
            m (int, default=None): If nmodes>1, this is the mode number to take 
                the gradient of. If left to default, it will take the gradient
                of the sum of all the modes. This arg doesn't do anything if
                nmodes=1.
        Returns:
            grad_f (dict of np.ndarray's): depends on geometry
                if cylindrical: {'R': R component, 'Z': Z component}
                if spherical: {'R': R component, 'Theta': Theta component}
                if cartesian: {'X': X component, 'Z': Z component}
        '''
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return
        
        if self._keff.size == 1 or m is None:
            ans = field_calcs.grad(self.field, self._coords[0], self._coords[1], geom=self._geom)
        else:
            m = int(m)
            ans = field_calcs.grad(self.modes[:,:,m], self._coords[0], self._coords[1], geom=self._geom)
        return ans

    def laplacian(self, m=None):
        ''' laplacian of field in poloidal plane

        This method takes the laplacian of the field in the correct geometry. Currently,
        this method is disabled for cases with phi as a np.ndarray.

        Args:
            m (int, default=None): If nmodes>1, this is th emode number to take
                the gradient of. If left to default, it will take the gradient
                of the sum of all modes. This arg doesn't do anything if nmodes=1.

        Returns:
            lap_f (np.ndarray): laplacian of field or mode component
        '''
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
            ans = field_calcs.laplacian(self.field, self._coords[0], self._coords[1], geom=self._geom)
        else:
            m = int(m)
            ans = field_calcs.laplacian(self.modes[:,:,m], self._coords[0], self._coords[1], geom=self._geom)
        return ans

    def plot(self, m=None, title=None, block=True):
        ''' quick and handy plotting function

        This method provides a quick filled contour plot of the field. Currently, this
        method is disabled for cases with phi as a np.ndarray.

        Args:
            m (int, default=None): If nmodes>1, this is the mode number to plot. If left
                to default, it will plot the sum of all modes. This arg doesn't do anything
            if nmodes=1.
            title (str, default=None): a title, if you want
            block (bool, default=True): if true, the matplotlib figure will block further
                execution. if false, it won't but you'll need a block eventually
        ''' 
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
            f, ax = plotting.contourf(self.field, self._grid, title=title)
        else:
            m = int(m)
            if title is None:
                title = 'm={0}'.format(self._keff[m])
            else:
                title += ' m={0}'.format(self._keff[m])
            f, ax = plotting.contourf(self.modes[:,:,m], self._grid, title=title)
        
        plt.show(block=block)

class VectorField(object):
    ''' VectorField class for manipulating nimpy scalar fields

    This class allows for manipulation of vector fields created by nimpy. It
    is an attribute of the NimData class for all instances of vector fields.
    Each component of the vector field is itself an instance of ScalarField.

    Attributes:
        Depends on geometry:
            cylindrical-->
                R (ScalarField): R component
                Phi (ScalarField): Phi component
                Z (ScalarField): Z component
            spherical-->
                R (ScalarField): R component
                Phi (ScalarField): Phi component
                Theta (ScalarField): Theta component
            cartesian-->
                X (ScalarField): X component
                Y (ScalarField): Y component
                Z (ScalarField): Z component
        Psi (ScalarField): streamfunction of the field. This is not included
            if phi is a np.ndarray (i.e. if the field is 3D)
    '''
    def __init__(self, f1, f2, f3, grid, geom, keff, phi):
        ''' VectorField initialization routine

        This method initializes an instance of VectorField. Based on
        the geometry, it assigns ScalarField instances for all three
        components of the field. Also, if the field is not 3D (i.e.
        phi is not a np.ndarray) it will calcuate poloidal streamlines.

        Args: f1 (np.ndarray): The R or X component of the field
              f2 (np.ndarray): The Phi or Y component of the field
              f3 (np.ndarray): The Z or Theta component of the field
              grid (tuple of np.ndarray's): the plotting grid for the data
              geom (str): geometry to use for calculations. Can be
                'cylindrical', 'spherical' or 'cartesian'
              keff (np.ndarray): the effective toroidal mode numbers
                (use [0] for axisymmetric case)
              phi (float or np.ndarray): the toroidal angle(s) to cast the
                field to if nmodes>1
        '''
        self._grid = grid
        self._geom = geom
        if self._geom == 'spherical':
            self._coords = grid_geom.convert_grid_to_spherical(grid[0], grid[1])
        else:
            self._coords = self._grid
        self._keff = keff
        self._phi = phi
        self._org_f1 = f1
        self._org_f2 = f2
        self._org_f3 = f3
        
        if self._geom == 'cylindrical':
            self._comps = ['R', 'Phi', 'Z']
        elif self._geom == 'spherical':
            self._comps = ['R', 'Phi', 'Theta']
        else:
            self._comps = ['X', 'Y', 'Z']

        fs = ['_f1','_f2','_f3']
        ms = ['_m1','_m2','_m3']
        for i, f in enumerate([f1, f2, f3]):
            setattr(self, self._comps[i], ScalarField(f, self._grid, self._geom, self._keff, self._phi))
            setattr(self, fs[i], getattr(self, self._comps[i]).field)
            if hasattr(getattr(self, self._comps[i]), 'modes'):
                setattr(self, ms[i], getattr(self, self._comps[i]).modes)

        self._x1, self._x3 = getattr(self, self._comps[0])._coords

        if type(self._phi) is not np.ndarray:
            self.Psi = ScalarField(self.__calc_streamlines__(), self._grid, self._geom, self._keff, self._phi)
    
    def __calc_streamlines__(self):
        if type(self._phi) is np.ndarray:
            print('method disabled for 3D fields. Recast to single toroidal angle and try again.')
            return
        else:
            return field_calcs.calc_poloidal_stream_func(self._f3, self._x1, x3=self._x3, geom=self._geom)
    def __change_phi__(self, new_phi):
        if hasattr(self, 'Psi'):
            del self.Psi
        self.__init__(self._org_f1, self._org_f2, self._org_f3, self._grid, self._geom, self._keff, new_phi)
    def __str__(self):
        line1 = 'VectorField nimpy object'
        if not any([np.any(self._f1), np.any(self._f2), np.any(self._f3)]):
            line2 = '\tempty array'
        else:
            line2 = ''
            for c in self._comps:
                min,max,avg,shape=getattr(self, c).__eval_field__()
                line2 += '{0}-comp:\n\t\tmin: {1}, max: {2}, avg: {3}, shape: {4}\n\t'.format(c, min, max, avg, shape)
            if hasattr(self, 'Psi'):
                min,max,avg,shape=self.Psi.__eval_field__()
                line2 += 'Psi:\n\t\tmin: {0}, max: {1}, avg: {2}, shape: {3}\n'.format(min, max, avg, shape)
        return line1+'\n\t'+line2

    def div(self, m=None):
        ''' divergence of field (no d/dphi or d/dy)

        This method takes the divergence of the field in the correct geometry. Currently,
        this method is disabled for cases with phi as a np.ndarray.

        Args:
            m (int, default=None): If nmodes>1, this is the mode number to take 
                the divergence of. If left to default, it will take the divergence
                of the sum of all the modes. This arg doesn't do anything if
                nmodes=1.
        Returns:
            div_f (np.ndarray): divergence of field
        '''
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
            ans = field_calcs.div(self._f1, self._f3, self._x1, self._x3, geom=self._geom)
        else:
            m = int(m)
            ans = field_calcs.div(self._m1[:,:,m], self._m3[:,:,m], self._x1, self._x3, geom=self._geom)
        return ans

    def curl(self, m=None):
         '''curl of field (no d/dphi or d/dy)

         This method takes the curl of the field in the correct geometry. Currently,
         this method is disabled for cases with phi as a np.ndarray.

         Args:
             m (int, default=None): If nmodes>1, this is the mode number to take 
                 the curl of. If left to default, it will take the curl
                 of the sum of all the modes. This arg doesn't do anything if
                 nmodes=1.
         Returns:
             curl_f (dict of np.ndarray's): curl of field
                 keys are geometry dependent:
                     'cylindrical'--> {'R', 'Phi', 'Z'}
                     'spherical'--> {'R', 'Phi', 'Theta'}
                     'cartesian'--> {'X', 'Y', 'Z'}
         '''
         if type(self._phi) is np.ndarray:
             print('method is disabled for 3D fields. Cast to single toroidal angle first')
             return

         if self._keff.size == 1 or m is None:
             ans = field_calcs.curl(self._f1, self._f2, self._f3, self._x1, self._x3, geom=self._geom)
         else:
             m = int(m)
             ans = field_calcs.curl(self._m1[:,:,m], self._m2[:,:,m], self._m3[:,:,m], self._x1, self._x3, geom=self._geom)
         return ans       

    def vec_laplacian(self, m=None):
         '''vector laplacian of field (no d/dphi or d/dy)

         This method takes the vector laplacian of the field in the correct geometry. Currently,
         this method is disabled for cases with phi as a np.ndarray.

         Args:
             m (int, default=None): If nmodes>1, this is the mode number to take 
                 the vector laplacian of. If left to default, it will take the vector
                 laplacian of the sum of all the modes. This arg doesn't do anything if
                 nmodes=1.
         Returns:
             vec_lap_f (dict of np.ndarray's): vector laplacian of field
                 keys are geometry dependent:
                     'cylindrical'--> {'R', 'Phi', 'Z'}
                     'spherical'--> {'R', 'Phi', 'Theta'}
                     'cartesian'--> {'X', 'Y', 'Z'}
         '''
         if type(self._phi) is np.ndarray:
             print('method is disabled for 3D fields. Cast to single toroidal angle first')
             return

         if self._keff.size == 1 or m is None:
             ans = field_calcs.vec_laplacian(self._f1, self._f2, self._f3, self._x1, self._x3, geom=self._geom)
         else:
             m = int(m)
             ans = field_calcs.vec_laplacian(self._m1[:,:,m], self._m2[:,:,m], self._m3[:,:,m], self._x1, self._x3, geom=self._geom)
         return ans 

    def f_dot_grad_f(self, m=None):
         '''convective derivative of field (no d/dphi or d/dy)

         This method takes the convective derivative of the field in the correct geometry. Currently,
         this method is disabled for cases with phi as a np.ndarray.

         Args:
             m (int, default=None): If nmodes>1, this is the mode number to take 
                 the convective derivative of. If left to default, it will take the convective
                 derivative of the sum of all the modes. This arg doesn't do anything if
                 nmodes=1.
         Returns:
             f_dot_grad_f (dict of np.ndarray's): convective derivative of field
                 keys are geometry dependent:
                     'cylindrical'--> {'R', 'Phi', 'Z'}
                     'spherical'--> {'R', 'Phi', 'Theta'}
                     'cartesian'--> {'X', 'Y', 'Z'}
         '''
         if type(self._phi) is np.ndarray:
             print('method is disabled for 3D fields. Cast to single toroidal angle first')
             return

         if self._keff.size == 1 or m is None:
             ans = field_calcs.f_dot_grad_f(self._f1, self._f2, self._f3, self._x1, self._x3, geom=self._geom)
         else:
             m = int(m)
             ans = field_calcs.f_dot_grad_f(self._m1[:,:,m], self._m2[:,:,m], self._m3[:,:,m], self._x1, self._x3, geom=self._geom)
         return ans 

    def plot(self, m=None, title=None, block=True):
        ''' quick and handy plotting function

        This method provides a quick filled contour plot of the field components. Currently, 
        this method is disabled for cases with phi as a np.ndarray.

        Args:
            m (int, default=None): If nmodes>1, this is the mode number to plot. If left
                to default, it will plot the sum of all modes. This arg doesn't do anything
            if nmodes=1.
            title (str, default=None): a title, if you want
            block (bool, default=True): if true, the matplotlib figure will block further
                execution. if false, it won't but you'll need a block eventually
        ''' 
        if type(self._phi) is np.ndarray:
            print('method is disabled for 3D fields. Cast to single toroidal angle first')
            return

        if self._keff.size == 1 or m is None:
            if not hasattr(self, 'Psi'):
                psi = None
            else:
                psi = self.Psi.field
            plotting.vec_contourf(self._f2, np.sqrt(self._f1**2+self._f3**2), self._grid, psi=psi, title=title)
        else:
            m = int(m)
            if title is None:
                title = 'm={0}'.format(self._keff[m])
            else:
                title += ' m={0}'.format(self._keff[m])
            plotting.vec_contourf(self._m2[:,:,m], np.sqrt(self._m1[:,:,m]**2+self._m3[:,:,m]**2), self._grid, title=title)
        plt.tight_layout() 
        plt.show(block=block)
