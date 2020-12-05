import numpy as np

'''
this module handles all the manipulation necessary for multiple-toroidal
mode nimrod runs. It depends on:
    numpy (as np)

In nimrod, fields with multiple toroidal modes are expressed as complex
with the following formula for total field:

    A(R,Z) = A(R,Z,m=0) + 2 * sum_(i=1->nmodes) Re[A(R,Z,m=i)]*sin(keff[i]*phi) - Im[A(R,Z,m=i)]*cos(keff[i]*phi)
'''
def complexify(in_dict):
    '''
    parses through a input dictionary from a dump .vtk file and casts toroidal
    direction as imaginary numbers in a single array.
    '''
    keys = list(set([k.replace('_2','') for k in list(in_dict.keys()) if '_2' in k]))
    out_dict = {}
    for k in keys:
        if '_' not in k:
            out_dict[k] = in_dict[k]+1j*in_dict[k+'_2']
        else:
            label = k.split('_')[0]+'_2_'+k.split('_')[1]
            out_dict[k] = in_dict[k]+1j*in_dict[label]
    leftovers = [l for l in list(in_dict.keys()) if '_2' not in l and l not in keys]
    for lf in leftovers:
        out_dict[lf] = in_dict[lf]

    return out_dict, keys

def combine_mode_dicts(big_list, m_list, cmpx_keys):
    '''
    input: big_list (list of dicts) a list of .vtk file dicts
           m_list (list of ints) a list of m values
           cmpx_keys (list of str) a list of keys of complex fields in the dict

    combines multiple mode dictionaries into one
    '''
    out_dict = {}
    sh = np.array(m_list).shape 
    for k in list(big_list[0].keys()):
        if k in cmpx_keys:
            out_dict[k] = np.zeros((big_list[0][k].shape + sh), dtype=np.complex128)
        else:
            out_dict[k] = big_list[0][k]
    
    for i,m in enumerate(m_list):
        for k in cmpx_keys:
            if len(out_dict[k].shape) is 3:
                out_dict[k][:,:,i] = big_list[i][k]
            else:
                out_dict[k][:,i] = big_list[i][k]
    
    return out_dict

def eval_complex_field(f, k_eff, phi=0, sum=False):
    '''
    Evaluate a complex field using the NIMROD convention:
    A(R,Z,Phi) = A_m=0 + 2*Sum(k=1,2,3...)[Re(A_m=k)*cos(k_eff*phi) - Im(A_m=k)*sin(k_eff*phi)
    where k_eff = m if geom='tor' and k_eff = m/2Pi if geom='rect'
    
    input: f (np.ndarray dim=2) input complex field
           k_eff (int, float, list) k_eff to be used
           phi (int, float, np.ndarray dim=1) the toroidal angle(s) to cast field to
           sum (bool, default=False) if true, all the modes will be summed, else the array will be returned with axis=2 
                representing the different modes

    output: f_out (np.ndarray dim=3 or 4) dim=3 if phi=float or int, dim=4 if phi=np.ndarray
                if sum is True, dim -= 1
    '''
    if f.dtype not in [complex, 'complex64', 'complex128']  or len(f.shape) != 3:
        print('Input field is not a complex field with multiple toroidal modes. Try again.')
        return f

    if type(phi) in [int, float, np.float64]:
        f_out = np.zeros(f.shape, dtype=float)
        for i, m in enumerate(k_eff):
            if m == 0.:
                f_out[:, :, i] = np.real(f[:, :, i])
            else:
                f_out[:, :, i] = 2.*np.real(f[:, :, i])*np.cos(m*phi) - 2.*np.imag(f[:, :, i])*np.sin(m*phi)
    else:
        f_out = np.zeros(f.shape+phi.shape, dtype=float)
        for i, m in enumerate(k_eff):
            if m == 0:
                f_out[:,:,i,:] = np.real(f[:,:,i])[:, :, np.newaxis]*np.ones_like(phi)
            else:
                f_out[:,:,i,:] = 2.*np.real(f[:, :, i])[:, :, np.newaxis]*np.cos(m*phi) - 2.*np.imag(f[:, :, i])[:, :, np.newaxis]*np.sin(m*phi)

    if sum:
        f_out = np.sum(f_out, axis=2)
        
    return f_out
