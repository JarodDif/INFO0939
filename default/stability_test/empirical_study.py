import subprocess
import struct
import numpy as np


rms_tab = []
dt_tab  = [*np.logspace(-8, -5, 3), *np.linspace(1e-5, 1e-4, 1000)]

for i, dt in zip(range(len(dt_tab)), dt_tab):

    with open("param_3d.txt", "r") as f:
        lines = f.readlines()

    lines[1] = str(dt) + "\n"
    lines[2] = str(5e2 * dt) + "\n"

    # Write the modified list back to the file
    with open("param_3d.txt", "w") as f:
        f.writelines(lines)

    result = subprocess.run(["../../fdtd", "param_3d.txt"], capture_output=True, text=True)

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
                rms = np.sqrt(np.sum(vals**2))/np.prod(vals.shape)
                rms_tab.append(rms)
                print(f"{dt = }\t: {rms = }")
                break
            
        if rms_tab[-1] > 1e3 :
            break
            
print(dt_tab[:len(rms_tab)], rms_tab)