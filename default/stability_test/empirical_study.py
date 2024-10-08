import subprocess
import struct
import numpy as np


h       = 10e-3
c       = 340
s_lim   = 3**(-1/2)
dt_lim  = s_lim * h / c
out_file = "empirical_study_"+str(h).replace(".", "_")

rms_tab = []
dt_tab  = np.array([*np.geomspace(dt_lim/1000, dt_lim/100, 2), *np.geomspace(dt_lim/10, dt_lim*10, 200)])


with open("param_3d.txt", "r") as f:
    lines = f.readlines()

lines[0] = f"{h:.2e}\n"

with open("param_3d.txt", "w") as f:
    f.writelines(lines)


for dt in dt_tab:

    with open("param_3d.txt", "r") as f:
        lines = f.readlines()
        
    lines[1] = f"{dt:.2e}\n"
    lines[2] = f"{500*dt:.2e}\n"
    lines[3] = "450\n"

    with open("param_3d.txt", "w") as f:
        f.writelines(lines)

    result = subprocess.run(["../fdtd", "param_3d.txt"], capture_output=True, text=True)

    print(result.stdout)

    def read_int(f):
        return struct.unpack("i", f.read(4))[0]

    def read_double(f):
        return struct.unpack("d", f.read(8))[0]

    with open("./out_p.dat", "rb") as f:
        
        nx = read_int(f)
        ny = read_int(f)
        nz = read_int(f)
        
        xmin = read_double(f)
        xmax = read_double(f)
        
        ymin = read_double(f)
        ymax = read_double(f)
        
        zmin = read_double(f)
        zmax = read_double(f)
        
        vals  = np.empty((nx, ny, nz))
            
        while True:
            
            try:
                index = read_int(f)
                time  = read_double(f)
                for _nz in range(nz):
                    for _ny in range(ny):
                        for _nx in range(nx):
                            vals[_nx, _ny, _nz] = read_double(f)
                
            except Exception:
                print(vals)
                rms_val = np.sqrt(np.sum(vals**2))/np.prod(vals.shape)
                rms_tab.append(rms_val)
                print(f"{dt = } and {dt_lim = }\t: {rms_val = }")
                break
            
        if rms_tab[-1] > 1e3 :
            break
            
np.savetxt(out_file, [dt_tab[:len(rms_tab)]*c/h, rms_tab])