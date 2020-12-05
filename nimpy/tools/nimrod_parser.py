import numpy as np
import struct
import subprocess
from os import devnull

'''
this module gets input and dump file information from a nimrod run
It depneds on:
    numpy (as np)
    struct
    subprocess
    os.devnull (as devnull)

The input parser looks at the nimset.out file created by nimrod which 
contains the entire namespace of user input for a nimrod run

The global dump file reader gets global values from an individual dumpfile
using binary reading operations

The dumpfile fields reader uses nimplot to create a vtk file of dump data
and then parses that as a text file to extract field information. This is
significantly less complicated than reading the fields straight from binary
and has the advantage of using nimplot which is kept up-to-date with nimrod
versions. 
'''
#--------------------------Input reader and functions----------------------#
def parse_line(a):
    ''' parses one line in nimset.out
    '''
    name = a.split('=')[0].strip().lower()
    alist = a.split('=')[1].split(',')
    if len(alist) > 2 or any([True for x in alist if '*' in x]):
        # This means it's an array #
        aout = []
        for x in alist:
            if '*' in x:
                for i in range(int(x.split('*')[0])):
                    try:
                        aout.append(float(x.split('*')[1]))
                    except ValueError:
                        aout.append(x.split('*')[1].replace('"','').strip())
            else: 
                try:
                    aout.append(float(x))
                except ValueError:
                    if '/n' not in x and x is not '':
                        aout.append(x.replace('"','').strip())
        if aout[-1] == '':
            aout = aout[0:-1]
        aout = np.array(aout)
    else:
        alist = alist[0]
        try:
            aout = int(alist)
        except ValueError:
            try:
                aout = float(alist)
            except ValueError:
                aout = alist.replace('"','')
                if aout == 'T':
                    aout = True
                elif aout == 'F':
                    aout = False
                else:
                    aout = aout
    return name, aout

def parse_nimset(fname='nimset.out'):
    ''' parses nimset.out created by nimrod
    
    this function creates a nested dictionary of the entire nimrod
    user input namespace

    Args:
        fname (str, default='nimset.out'): file to parse

    Returns:
        inputs (dict): a nested dictionary with sub dictionaries
            named after the input sections. All information on
            input values can be found in nimset/input.f
    '''
    dict_name_list = []
    input_dict = {}
    add_to = False
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            #print i+1
            if '&' in line:
                if line.split('&')[1].strip().lower() not in dict_name_list:
                    dict_name_list.append(line.split('&')[1].strip().lower())
                    add_to = True
                    a = {}
            if line.replace('/n','').strip() == '/':
                if add_to == True:
                    input_dict[dict_name_list[-1]] = a
                    del a
                    add_to = False
                else:
                    add_to = False
            if add_to and '&' not in line:
                name, aout = parse_line(line)
                a[name] = aout
    f.close()
    return input_dict

#-------------------------------Global dump file reader-----------------------#
def get_endian(fname):
    ''' gets the endian of a binary file
    '''
    fid = open(fname, 'rb')
    if struct.unpack('<' + 'i', fid.read(4))[0] != 8:
        endian = '>'
    else:
        endian = '<'
    fid.close()
    return endian

def read_header(fid, endian, intsize=4, output=False):
    ''' reads the header of a binary file
    '''
    val = int(struct.unpack(endian+'i',fid.read(intsize))[0])
    if output:
        return val
    else:
        return None

def read_value(fid, endian):
    ''' reads a value from a binary file
    '''
    read_header(fid, endian)
    val = struct.unpack(endian+'d',fid.read(8))[0]
    read_header(fid, endian)
    return val

def read_array(fid, endian, shape, floatsize=8):
    ''' reads an array from a binary file
    '''
    read_header(fid, endian)
    if (type(shape) == type(1)) or (type(shape) == type(1.)):
        vals = np.array(struct.unpack(endian+'d'*int(shape),fid.read(floatsize*shape)), dtype=float)
    else:
        shape = list(shape)
        shape.reverse()
        n_elements = 1
        for n_dim in shape:
            n_elements *= n_dim
        vals = np.array(struct.unpack(endian+'d'*int(n_elements),fid.read(floatsize*n_elements)), dtype=float).reshape(shape).T
    read_header(fid, endian)
    return vals

def read_global_vals(fname):
    ''' reads in global values from nimrod dumpfile

    this function opens a nimrod dumpfile and reads
    it to get global values (not fields)

    Args:
        fname (str): dumpfile name

    Returns:
        global_vals (dict): a dictionary of the 
            global vals. With the following:
        'time' (float): time of dumpfile
        'step' (int): time step of dumpfile
        'nbl' (int): number of blocks
        'nrbl' (int): number of r blocks
        'poly_deg' (int): poly degree used
        'nmodes' (int): number of toroidal modes
        'keff' (np.ndarray): toroidal mode numbers
        'nseams' (int): number of seams
    '''
    endian = get_endian(fname)
    fid = open(fname, 'rb')
    
    time = read_value(fid, endian)
    step = int(read_value(fid, endian))
    nbl = int(read_value(fid, endian))
    nrbl = int(read_value(fid, endian))
    poly_deg = int(read_value(fid, endian))
    nmodes = int(read_value(fid, endian))
    keff = read_array(fid, endian, nmodes).astype(int)
    nseams = nbl+1

    fid.close()
    
    return {'time': time, 'step': step, 'nbl': nbl, 'nrbl': nrbl, 'poly_deg': poly_deg, 'nmodes': nmodes, 'keff': keff, 'nseams': nseams}

#-----------------------------------Dump file fields reader----------------------#
def file_len(fname):
    ''' counts the number of lines in a file
    '''
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def read_vec(fname):
    '''
    Reads in a .vtk file produced by nimplot with option '5' (vector plot). Currently, works only with one mode.
    dependencies: numpy (import numpy as np)
    
    :param fname: (str) .vtk name to read
    :return: (dict) dictionary containing grid and fields 
    
    '''
    f = open(fname)
    vecname = []
    scalname = []
    for i in range(file_len(fname)):
        next = f.readline()
        if any([x.isalpha() for x in next]) and any([y for y in [x for x in next if x.isalpha()] if y!='E']):
            if 'POINTS' in next:
                n = int(''.join([x for x in next if x.isdigit()]))
                BigDict = {}
                R = np.zeros(n)
                Z = np.zeros(n)
                for i in range(n):
                    line = [x for x in f.readline().replace('\n', '').split(' ') if x!='']
                    R[i] = float(line[0])
                    Z[i] = float(line[1])
                BigDict['R'] = R
                BigDict['Z'] = Z
            if 'VECTORS' in next:
                vecname.append(next.split(' ')[1])
                R = np.zeros(n)
                Z = np.zeros(n)
                Phi = np.zeros(n)
                for i in range(n):
                    line = [x for x in f.readline().replace('\n', '').split(' ') if x!='']
                    R[i] = float(line[0])
                    Z[i] = float(line[1])
                    Phi[i] = float(line[2])
                BigDict[vecname[-1]+'_R'] = R
                BigDict[vecname[-1]+'_Z'] = Z
                BigDict[vecname[-1]+'_Phi'] = Phi
            if 'SCALARS' in next:
                scalname.append(next.split(' ')[1])
                f.readline()
                X = np.zeros(n)
                for i in range(n):
                    X[i] = float(f.readline().replace('\n', '').replace(' ',''))
                BigDict[scalname[-1]] = X
    return BigDict

def prepare_nimplot_vec(dump_name, vtk_name, m):
    '''
    Runs nimplot silently for '5' (vector plot) option.
    dependencies: subprocess (import subprocess)
    
    :param dump_name: (str) dump file to run nimplot with (dump.#####) 
    :param vtk_name: (str) fname for nimplot output, must end with .vtk
    :return: none
    
    '''
    nim = "nimplot.temp.in"
    cmd = "touch {0} ".format(nim)
    subprocess.Popen(cmd.split())
    subprocess.Popen('cp $NIMSRCDIR/draw/drawvec.in .', shell=True)
    f = open(nim, 'w')
    f.write('blah\n')
    f.write('\n')
    f.write('5\n')
    f.write(dump_name + '\n')
    f.write('\n')
    f.write('o\n')
    if m is None or int(m) is 0:
        f.write('0\n')
        f.write('y\n')
    else:
        f.write('{0}\n'.format(int(m)))
    f.write('\n')
    f.write(vtk_name + '\n')
    f.write('\n')
    f.write('\n')
    f.write('0\n')
    f.close()
    print("Running nimplot, might take a while...")
    subprocess.Popen("nimplot < {0}".format(nim), shell=True, stdout=open(devnull, 'wb'),
                     stderr=subprocess.STDOUT).wait()
    subprocess.Popen("rm {0} drawvec.in".format(nim).split())


def parse_dump_fields(dump_name, vtk_name, vtk_write, m):
    '''
    Main parsing function for nimplot vector option.
    
    :param dump_name: (str) dump file to parse (dump.#####)
    :param vtk_name: (str) vtk name for nimplot to make (must end with .vtk)
    :param vtk_write: (bool) if True, will not delete vtk file after parsing
    :return: big_dict (dict) dictionary of parsed, interpolated grid and fields
    
    '''
    prepare_nimplot_vec(dump_name, vtk_name, m)
    big_dict = read_vec(vtk_name)
    if not vtk_write:
        subprocess.Popen("rm {0}".format(vtk_name).split())
    return big_dict
