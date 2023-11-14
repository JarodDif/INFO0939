import struct
import numpy as np

def read_int(f):
    return struct.unpack("i", f.read(4))[0]

def read_double(f):
    return struct.unpack("d", f.read(8))[0]

with open("./default/example_inputs/stability_test/out_p.dat", "rb") as f:
    
    nx = read_int(f)
    ny = read_int(f)
    nz = read_int(f)
       
    xmin = read_double(f)
    xmax = read_double(f)
    
    ymin = read_double(f)
    ymax = read_double(f)
    
    zmin = read_double(f)
    zmax = read_double(f)
    
    
    while True:
        
        vals  = np.empty((nx, ny, nz))
        try:
            index = read_int(f)
            time  = read_double(f)
            for _nz in range(nz):
                for _ny in range(ny):
                    for _nx in range(nx):
                        vals[_nx, _ny, _nz] = read_double(f)
            
        except Exception:
            print(np.sum(np.abs(vals))/np.prod(vals.shape))
            exit(-1)