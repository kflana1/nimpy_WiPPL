import numpy as np
import struct
import argparse
import h5py
import matplotlib.pyplot as plt
import os.path

'''
this module reads the energy.bin files created by nimrod and 
saves the data into a python-friendly hdf5 file.
It depends on:
    numpy (as np)
    struct
    argparse
    h5py
    matplotlib.pyplot (for plotting the result)
    os.path

The binary reading has been mostly taken from other nimrod
users.
'''
def read_energy(filename='energy.bin'):
    ''' reads energy.bin file to dictionary

    Using struct, this reads a energy.bin file
    and returns a dictionary of contained data

    Args:
        filename (str, default='energy.bin'):
            the energy.bin filename to read
    
    Returns:
        energy (dict): dictionary with keys:
            't' (np.ndarray): time array
            'kmode' (np.ndarray): toroidal mode array
            'magen' (np.ndarray): magentic energy
            'kinen' (np.ndarray): kinetic energy
    '''
    #0	istep
    #1	t
    #2	imode
    #3	k
    #4	E_m
    #5	E_k
    int_size =4
    float_size = 4
    #open file
    #first arg is filename, second is mode 'rb=read binary'
    #f=open('power_flux.bin.NOEQ','rb')
    f=open(filename,'rb')
    #go to the end of the file and find the size
    # seek: first arg is offset, second is to where: 0=begin, 1=current, 2=end
    f.seek(0,2)
    file_size=f.tell()
    #presuming that size is the number of data points in file
    #go to front of the file and find the number of lines per block
    f.seek(0)
    #first byte is length of current line
    #last byte is also length of current line
    temp=f.read(4)
    #llb=line length bytes
    blah=struct.unpack(">l",temp)
    llb=blah[0]
    #note line length is in bytes, divide by float_size to get 
    #number of data points per line
    #dppl=data points per line
    dppl=np.fix(llb/float_size)
    #declare workarray linein
    #array that stores one line of data read from file
    linein=np.zeros(int(dppl),dtype='f',order = 'F')
    #count the number of lines in a data block
    #lpb=lines per block
    lpb=0
    while blah[0] > 0:
    # maybe do a loop from i=0 to line_length, read in temp
    # then copy all to tf array? would be more generic I think
    #    for i in (0,dppl-1):
        temp=f.read(llb)
    #    tf = struct.unpack(">ffffff", temp)
        lpb=lpb+1
        temp=f.read(4)
        blah=struct.unpack(">l",temp)
        temp=f.read(4)
        blah=struct.unpack(">l",temp)
    #now figure out the number of blocks of data in the file
    #  note this assumes 4 byte data for all elements
    #  (dppl+2) includes the beginning and ending
    #  (dppl+2)*lpb incorporates all lines in the block
    #  (dppl+2)*lpb+2 incorporates the 0 0 at the end of a block 
    #dbs=data blocks
    dbs=file_size/(float_size*((dppl+2)*lpb+2))
    #print(dbs)
    #BEGIN READING IN ALL DATA HERE
    #declare variables
    time  = np.zeros(int(dbs), dtype = 'f',order = 'F')
    kmode = np.zeros(lpb, dtype = 'f',order = 'F')
    magen = np.zeros([lpb, int(dbs)], dtype = 'f',order = 'F')
    kinen = np.zeros([lpb, int(dbs)], dtype = 'f',order = 'F')
    f.seek(0)
    iy=0
    while iy < dbs:
        ix=0
        while ix < lpb:
            temp=f.read(4)
#        blah=struct.unpack(">l",temp)
            temp=f.read(6*float_size)
            linein = struct.unpack(">ffffff", temp)
            if ix==0:
                time[iy] = linein[1]
            if iy==0:
                kmode[ix] = linein[2]
            magen[ix,iy] = linein[4] 
            kinen[ix,iy] = linein[5] 
            temp=f.read(4)
#        blah=struct.unpack(">l",temp) 
            ix=ix+1
        temp=f.read(4)
#on last pass, this next one should not read in anything, but should be ok
        temp=f.read(4)
        iy=iy+1

    f.close()

    return {'t':time, 'kmode':kmode, 'magen':magen, 'kinen':kinen}

def combine(fnames=['energy.bin1', 'energy.bin']):
    ''' combines energy.bin information from multiple bin files

    restarting nimrod will make a new energy.bin file, so 
    this function will combine the data from multiple files into 
    one dictionary

    Args:
        fnames (list of str): list of energy.bin files
         NOTE: should be in chronological order (first is earliest)
    
    Returns:
        energy (dict): dictionary with keys:
            't' (np.ndarray): time array
            'kmode' (np.ndarray): toroidal mode array
            'magen' (np.ndarray): magentic energy
            'kinen' (np.ndarray): kinetic energy
    '''
    a = read_energy(filename=fnames[0])
    time = a['t']
    kmode = a['kmode']
    magen = a['magen']
    kinen = a['kinen']
    
    if len(fnames) == 1:
        return {'t':time, 'kmode':kmode, 'magen':magen, 'kinen':kinen}
    
    for fname in fnames[1::]:
        a = read_energy(fname)
        time = np.concatenate((time,a['t']))
        magen = np.concatenate((magen,a['magen']),axis=1)
        kinen = np.concatenate((kinen,a['kinen']),axis=1)

    return {'t':time, 'kmode':kmode, 'magen':magen, 'kinen':kinen}

def plot_en(fname='energy.bin.h5',saveit=None,block=True):
    f = h5py.File(fname,'r')
    time = f['t'].value*1.e3
    kmode = f['kmode'].value
    magen = f['magen'].value
    kinen = f['kinen'].value
    f.close()

    f,ax = plt.subplots(2,1,sharex=True,figsize=(10,9))
    ktot = np.zeros_like(time)
    mtot = np.zeros_like(time)
    for i,k in enumerate(kmode):
        ax[0].plot(time,magen[i,:],lw=3,label='m={0}'.format(int(k-1)))
        ax[1].plot(time,kinen[i,:],lw=3,label='m={0}'.format(int(k-1)))
        mtot += magen[i,:]
        ktot += kinen[i,:]
    ax[0].plot(time,mtot,lw=3,label='total')
    ax[1].plot(time,ktot,lw=3,label='total')
    ax[0].legend(loc='best', fontsize=22)
    ax[1].legend(loc='best', fontsize=22)
    ax[1].set_xlabel('t (ms)', fontsize=22)
    ax[0].set_title('Magnetic Energy', fontsize=24)
    ax[1].set_title('Kinetic Energy', fontsize=24)
    ax[0].set_xlim([0,time[-1]])
    ax[1].set_xlim([0,time[-1]])

    for axx in ax:
        for loc in ['top', 'right', 'bottom', 'left']:
            axx.spines[loc].set_linewidth(2)
        for tick in axx.xaxis.get_major_ticks():
            tick.label1.set_fontsize(16)
        for tick in axx.yaxis.get_major_ticks():
            tick.label1.set_fontsize(16)
        axx.xaxis.set_tick_params(width=3, length=8)
        axx.yaxis.set_tick_params(width=3, length=8)
                          
    plt.tight_layout()
    if saveit is None:
        plt.show(block=block)
    else:
        plt.savefig(saveit,bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Python program for reading the energy.bin time history created by nimrod.')
    parser.add_argument('-n', '--name', nargs='+', default=['energy.bin'], help='name of the energy.bin file[s] if you have more than one list them here in order from earliest time to latest')
    parser.add_argument('-p', '--plot', action='store_true', help='flag to plot the energy')
    args = parser.parse_args()
    energy = combine(args.name)
    if args.plot and os.path.isfile('energy.bin.h5'):
        plot_en()
    else:
        with h5py.File('energy.bin.h5', 'w') as f:
            f.create_dataset('t', data=energy['t'], compression='lzf')
            f.create_dataset('kmode', data=energy['kmode'])
            f.create_dataset('magen', data=energy['magen'], compression='lzf')
            f.create_dataset('kinen', data=energy['kinen'], compression='lzf')
        print('Energy data saved to energy.bin.h5')
        if args.plot:
            plot_en()
