import numpy as np
import struct

dims = (101, 2, 101)
lims = (0,1), (0,1), (0,1)

def write_int(f, n: int):
    f.write(n.to_bytes(4, byteorder="little"))

def write_double(f, d:float):
    f.write(struct.pack('d', d))

def write_dims(f, dims:tuple[int,int,int]):
    for d in dims:
        write_int(f, d)

def write_lims(f, lims:tuple[tuple[float, float], tuple[float, float], tuple[float, float]]):
    for (l,u) in lims:
        write_double(f, l)
        write_double(f, u)

# Will write for stepindex 0, timestep 0
def write_data(f, vals):
    write_int(f, 0)
    write_double(f, 0)
    for p in range(vals.shape[2]):
        for n in range(vals.shape[1]):
            for m in range(vals.shape[0]):
                write_double(f, vals[m,n,p])

vals = np.full(dims, 340)
vals[:, :, dims[2]//4-1:dims[2]//4+2]  = 0
vals[45, :, dims[2]//4-1:dims[2]//4+2] = 340
vals[55, :, dims[2]//4-1:dims[2]//4+2] = 340

with open("../example_inputs/young/in_c_young_abs.dat", "wb") as f:
    write_dims(f,dims)
    write_lims(f,lims)
    write_data(f,vals)
    
vals = np.full(dims, 340)
vals[:, :, dims[2]//4-1:dims[2]//4+2]  = -340
vals[45, :, dims[2]//4-1:dims[2]//4+2] = 340
vals[55, :, dims[2]//4-1:dims[2]//4+2] = 340

with open("../example_inputs/young/in_c_young_ref.dat", "wb") as f:
    write_dims(f,dims)
    write_lims(f,lims)
    write_data(f,vals)