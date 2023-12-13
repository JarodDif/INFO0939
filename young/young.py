import numpy as np
import struct
import matplotlib.pyplot as plt

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

vals                    = np.full(dims, 340)
vals[:, :, dims[2]//2]  = 0
vals[45, :, dims[2]//2] = 340
vals[55, :, dims[2]//2] = 340

im = plt.imshow(vals[:, :, dims[2]//2], aspect='auto')
plt.colorbar(im)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

indices = np.where(vals == 0)
print(indices)
ax.scatter(*indices, c="k")
plt.show()

with open("../example_inputs/simple3d/in_c_custom_young.dat", "wb") as f:
    write_dims(f,dims)
    write_lims(f,lims)
    write_data(f,vals)