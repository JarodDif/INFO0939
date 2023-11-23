import struct
import sys
import glob
import math
import numpy as np

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

if(len(sys.argv) < 2):
    print(f"USAGE : python {sys.argv[0]} parameter_file NX NY NZ")
    sys.exit()

param_file = sys.argv[1]



#read paramfile
with open(param_file, "r") as f:
    dx = float(f.readline())
    dt = float(f.readline())
    maxt = float(f.readline())
    outputrate = int(f.readline())
    [f.readline() for i in range(3)]
    first_output = f.readline()
    filename = first_output.split(" ")[2]
    filename = filename.strip("\r\n")

#find out how many timesteps there will be
timesteps = int(maxt / dt / outputrate) + 1

#find all corresponding files
files = glob.glob(f"*_{filename}")
file_handles = [None]*len(files)

#find back dimensions
dims = [0]*3

for fname in files:
    coords = [int(c) for c in fname.split("_")[:3]]
    for i in range(3):
        dims[i] = max(dims[i], coords[i])
dims = tuple(d+1 for d in dims)

sizes = np.zeros(dims+(3,), dtype=int)
dom_lims = np.zeros(dims+(3,2))

for i, fname in enumerate(files):
    coords = tuple(int(c) for c in fname.split("_")[:3])
    file_handles[i] = open(fname, "rb")
    sizes[coords] = [read_int(file_handles[i]) for j in range(3)]
    dom_lims[coords] = [[read_double(file_handles[i]), read_double(file_handles[i])] for j in range(3)]

NX = sizes[:,0,0,0].sum()
NY = sizes[0,:,0,1].sum()
NZ = sizes[0,0,:,2].sum()
N = (NX, NY, NZ)

xmin, xmax = dom_lims[:,0,0,0].flatten()[[0,-1]]
ymin, ymax = dom_lims[0,:,0,1].flatten()[[0,-1]]
zmin, zmax = dom_lims[0,0,:,2].flatten()[[0,-1]]

domain = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]

#compute limits
lims = [x[:] for x in [[None]*3]*len(files)]
lims = np.zeros((len(files), 3, 3), dtype=int)
for i, fname in enumerate(files):
    coords = [int(c) for c in fname.split("_")[:3]]
    for j in range(3):
        s = N[j]*coords[j]//dims[j]
        e = N[j]*(coords[j]+1)//dims[j] -1
        n = e-s+1
        lims[i,j] = s,e,n    

print("Grid :", N)
print("Domain : ", domain)

output_file = open(filename, "wb")

for n in N:
    write_int(output_file, int(n))
for (_m, _M) in domain:
    write_double(output_file, _m)
    write_double(output_file, _M)

#prepare data for one time step
data = np.zeros((NX,NY,NZ))

for t in range(timesteps):
    tstep, time = 0, 0.0
    for handle, lim in zip(file_handles, lims):
        tstep, time = read_int(handle), read_double(handle)
        (sx, ex, nx), (sy, ey, ny), (sz, ez, nz) = lim
        for p in range(nz):
            for n in range(ny):
                for m in range(nx):
                    data[p+sz, n+sy, m+sx] = read_double(handle)
    write_int(output_file, tstep)
    write_double(output_file, time)
    for p in range(NZ):
        for n in range(NY):
            for m in range(NZ):
                write_double(output_file, float(data[p,n,m]))
    print(f"{t+1:3d} / {timesteps:3d}")

output_file.close()
for h in file_handles:
    h.close()