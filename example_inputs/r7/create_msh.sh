ml Info0939Tools

data2gmsh -cutz 0.2 0_0_0_out_p.dat
mv cutz.msh 0_0_0_cutz.msh
data2gmsh -cutz 0.2 0_1_0_out_p.dat
mv cutz.msh 0_1_0_cutz.msh
data2gmsh -cutz 0.2 1_0_0_out_p.dat
mv cutz.msh 1_0_0_cutz.msh
data2gmsh -cutz 0.2 1_1_0_out_p.dat
mv cutz.msh 1_1_0_cutz.msh

rm *_out_p.dat