import numpy as np
import argparse
import subprocess
from os import path
import warnings
from nimpy.tools.nimrod_parser import parse_nimset, read_global_vals, parse_dump_fields
from nimpy.tools.grid_geom import make_regular_grid
from nimpy.tools.toroidal_modes import complexify, combine_mode_dicts
from nimpy.tools.fileIO import save_nimpy_dict

'''
nimh5 is used to create hdf5 files from nimrod binary dump files.
It should be used as a function, not a module.
Try: nimh5 --help
'''
def read_single_mode(fname, m=0, n=None):
    big_dict = parse_dump_fields(fname, fname+'.vtk', False, m)
    big_dict = make_regular_grid(big_dict, n=n)
    return big_dict

def proc_one_dump(fname, n, m_list=None):

    big_dict = {}
    global_vals = read_global_vals(fname)
    if m_list is None:
        m_list = list(range(int(global_vals['nmodes'])))

    if m_list == [0]:
        a_dict = read_single_mode(fname, m=0, n=n)
        for k in [x for x in list(a_dict.keys()) if '_2' not in x]:
            big_dict[k] = a_dict[k]
    else:
        big_list = []
        for i, m in enumerate(m_list):
            print(("processing toroidal mode m={0}".format(m)))
            a_dict, cmpx_keys = complexify(read_single_mode(fname, m=m, n=n))
            big_list.append(a_dict)
        big_dict = combine_mode_dicts(big_list, m_list, cmpx_keys)

    big_dict['global_vals'] = global_vals

    return big_dict

def parse_dump_list(dump_list=None, dump_all=False):

    if dump_all:
        out_list = subprocess.Popen("ls -1v dump.* | grep -v *.h5", shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0:-1]
    else:
        if dump_list is None:
            out_list = subprocess.Popen("ls -1v dump.* | grep -v *.h5 | tail -n 1", shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0:-1]
        else:
            out_list = []
            for d in dump_list:
                if 'dump.' not in d:
                    d = "dump.{0:05d}".format(int(d))
                if path.isfile(d):
                    out_list.append(d)
            if len(out_list) != len(dump_list):
                warnings.warn('Not all requested dumps are in this directory. Please check and try again')
    return out_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Python program for running nimplot and producing python-readable datafiles. Filetype is hdf5 (.h5) which can easily be read in as dictonaries in python using h5py package.')
    parser.add_argument('-d', '--dump_list', nargs='+', default=None, help='Dump file(s) to process. Can input full file name or just number. If not called, will automatically use dump with largest number.')
    parser.add_argument('-a', '--dump_all', action='store_true', help='Flag to process all dump files in this directory. This takes a while, so you might want to use nohup or screen.')
    parser.add_argument('-n', '--n_points', nargs='+', default=None, type=int, help='Number of grid points to use when interpolating onto regular grid. If one number is given, the grid will be square. If two numbers given, grid will be a rectangle. If not called, the default will be to return a square regular grid with approximately as many points as the nimrod mesh used.')
    parser.add_argument('-m', '--m_list', nargs='+', default=None, type=int, help='Toroidal mode number(s) to process for each dump file. If numbers are entered, then only those modes will be processed. If not called, all available modes in the dump files will be processed.') 
    parser.add_argument('-i', '--input_vals', action='store_true', help='Flag to store a dictionary of input vals in hdf5 file. This will parse nimset.out for values.')
    args = parser.parse_args()
    dump_list = parse_dump_list(dump_list=args.dump_list, dump_all=args.dump_all)
    if args.input_vals:
        if path.isfile('nimset.out'):
            input_vals = parse_nimset()
        else:
            raise IOError('nimset.out not found')
    for d in dump_list:
        print(('{0}'.format(d)))
        big_dict = proc_one_dump(d, args.n_points, args.m_list)
        if args.input_vals:
            big_dict['input_vals'] = input_vals
        save_nimpy_dict(d+'.h5', big_dict)
        print('--------------------------------')
   
