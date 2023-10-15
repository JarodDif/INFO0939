import struct

import numpy as np
import matplotlib.pyplot as plt

def read_int(f):
    b = f.read(4)
    if len(b) != 4:
        raise Exception("Not enough bytes")
    return int.from_bytes(b, byteorder="little")

def read_double(f):
    b = f.read(8)
    if len(b) != 8:
        raise Exception("Not enough bytes")
    return struct.unpack("d", b)[0]


with open("mpi/example_inputs/simple3d/in_c_custom.dat", "rb") as f:
    dims = read_int(f), read_int(f), read_int(f)
    print("Dimensions :", dims)

    lims = (read_double(f), read_double(f)), (read_double(f), read_double(f)), (read_double(f), read_double(f))
    print("Limits :", lims)

    #read each time step
    while True:
        stepindex = -1
        try:
            stepindex = read_int(f)
        except:
            break
        timestep = read_double(f)
        print(f"Step {stepindex} at t={timestep}")
        vals = np.zeros(dims)
        for p in range(dims[2]):
            for n in range(dims[1]):
                for m in range(dims[0]):
                    vals[m,n,p] = read_double(f)
        print(vals)

        ax = plt.figure().add_subplot(projection="3d")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        for i, v in np.ndenumerate(vals):
            coords = [
                l + n/(d-1) * (u-l) for n,d,(l,u) in zip(i,dims,lims)
            ]
            ax.text(*coords, str(v))

        plt.show()
