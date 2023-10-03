import numpy as np
import struct

dims = (2,2,2)
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

vals = np.zeros(dims)

for n in range(dims[1]):
    for p in range(dims[2]):
        vals[0,n,p] = 170
        vals[1,n,p] = 340

with open("example_inputs/simple3d/in_c_custom.dat", "wb") as f:
    write_dims(f,dims)
    write_lims(f,lims)
    write_data(f,vals)