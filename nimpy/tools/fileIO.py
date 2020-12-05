import numpy as np
import h5py
import sys
from io import StringIO

'''
nimpy.file_io is the nimpy module for preparing the hdf5/dictionary structure
this module depends on:
    numpy (defined as np)
    h5py
    sys
    cStringIO (defined as StringIO)
'''
def recursive_save_dict_to_h5(h5, path, dic):
    ''' function used in save_dict_to_h5 in order to get recursion
    '''
    for key, item in list(dic.items()):
        if type(item) in [np.ndarray, np.generic]:
            h5.create_dataset(path+key, data=item, compression='lzf')
        elif type(item) != dict:
            try:
                h5.create_dataset(path+key, data=item)
            except TypeError:
                raise ValueError('Cannot save %s type'%type(item))
        else:
            recursive_save_dict_to_h5(h5, path+key+'/', item)

def save_dict_to_h5(fname, dic):
    '''save a dictionary of dictionaries (of arbitrary depth)

    saves a dictionary of arbitrary depth to an hdf5 file
    using the h5py package. Makes calls to 
    recursive_save_dict_to_h5

    Args:
        fname (str): name of h5 file you want to make
        dic (dict): dictionary you want to save
    '''
    with h5py.File(fname, 'w') as h5:
        recursive_save_dict_to_h5(h5, '/', dic)

def recursive_load_dict_from_h5(h5, path):
    ''' function used in load_h5_to_dict in order to get recursion
    '''
    out_dict = {}
    for key, item in list(h5[path].items()):
        if type(item) == h5py._hl.dataset.Dataset:
            out_dict[key] = item.value
        elif type(item) == h5py._hl.group.Group:
            out_dict[key] = recursive_load_dict_from_h5(h5, path+key+'/')
    return out_dict


def load_h5_to_dict(fname):
    '''load a dictionary from a hdf5 file

    loads a dictionary from an arbitrary hdf5 file
    using h5py package. Makes calls to 
    recursive_load_dict_from_h5

    Args:
        fname (str): filename of hdf5 file you
            wish to read in. (must include extension)
    '''
    with h5py.File(fname, 'r') as h5:
        return recursive_load_dict_from_h5(h5, '/')

def dic_length(dic):
    '''gets the number of values in recursive dictionary

    recursively counts all values in a dictionary that
    are not themselves dictionaries

    Args:
        dic (dict): dictionary you want to know the
            length of
    
    Returns:
        count (int): number of values in dic
    '''
    count = 0
    for k in list(dic.keys()):
        if type(dic[k]) is dict:
            count += dic_length(dic[k])
        else:
            count += 1
    return count

def dic_print(dic, indent=0):
    '''prints the values in a nested dictionary

    prints the values in a nested dictionary in a
    clear way

    Args:
        dic (dict): dictionary to print
        indent (int, default=0): initial
            tab indent of print
    '''
    for key, value in dic.items():
        print('\t' * indent + str(key))
        if type(value) is dict:
            dic_print(value, indent+1)
        else:
            if type(value) in [str, int, float, bool, list, np.int64, np.float64]:
                print('\t' * (indent+1) + '{0}'.format(value))
            elif type(value) in [np.ndarray, np.generic]:
                if value.size < 21 or type(value[0]) is str:
                    print('\t' * (indent+1) + '{0}'.format(value))
                elif np.any(value):
                    print('\t' * (indent+1) + 'size: {0}, min: {1}, max:{2}'.format(value.shape, np.nanmin(value), np.nanmax(value)))
                else:
                    print('\t' * (indent+1) + 'empty array')

def dic_print_to_str(dic, indent=0):
    '''returns a string printing nested dictionary conents

    calls dic_print, but instead of printing, returns a 
    string

    Args:
        dic (dict): dictionary to print to string
        indent (int, default=0): initial
            tab indent of print
    
    Returns:
        dic_string (str): string of dicitionary values
    '''
    old_stdout = sys.stdout
    sys.stdout = ans = StringIO()
    dic_print(dic, indent=indent)
    sys.stdout = old_stdout
    return ans.getvalue()

def create_nimpy_dict(d):
    ''' creates the nested nimpy dictionary from a flat dictionary

    this function organizes the nimpy dictionary structure that is
    used throughout nimpy

    Args: 
        d (dict): a flat (i.e. not nested) dictionary from output
            of parsing a dump file
    
    Returns:
        nimpy_dict (dict): nested dictionary with nimpy data
            structure
    '''
    ans = {}
    
    ans['grid'] = {'R': d['R'], 'Z': d['Z']}

    vec_list = []
    scal_list = []
    for k in list(d.keys()):
        if type(d[k]) is dict:
            ans[k] = d[k]
        elif '_R' in k or '_Z' in k or '_Phi' in k:
            name = k.replace('_R','').replace('_Phi','').replace('_Z','')
            vec_list.append(name)
        elif k not in ['R', 'Z']:
            scal_list.append(k)
    vec_list = list(set(vec_list))

    ans['eq_fields'] = {'vector':{}, 'scalar':{}}
    ans['fields'] = {'vector':{}, 'scalar':{}}

    for v in vec_list:
        if '0' in v:
            ans['eq_fields']['vector'][v] = {'R': d[v+'_R'], 'Phi': d[v+'_Phi'], 'Z': d[v+'_Z']}
        else:
            ans['fields']['vector'][v] = {'R': d[v+'_R'], 'Phi': d[v+'_Phi'], 'Z': d[v+'_Z']}

    for s in scal_list:
        if '0' in s:
            ans['eq_fields']['scalar'][s] = d[s]
        else:
            ans['fields']['scalar'][s] = d[s]
    
    return ans

def save_nimpy_dict(fname, d):
    ''' wrapper function for creating and saving a nimpy dictionary
    
    calls create_nimpy_dict and save_dict_to_h5 on a flat (i.e. not
    nested) dictionary

    Args:
        fname (str): filename to write dictionary to
        d (dict): a flat (i.e. not nested) dictionary from output
            of parsing a dump file
    '''
    a = create_nimpy_dict(d)
    save_dict_to_h5(fname, a)
