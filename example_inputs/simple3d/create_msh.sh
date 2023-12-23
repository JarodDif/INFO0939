ml Info0939Tools

data2gmsh -cuty 0.5 0_0_0_out_p.dat
mv cuty.msh 0_0_0_cuty.msh
data2gmsh -cuty 0.5 1_0_0_out_p.dat
mv cuty.msh 1_0_0_cuty.msh

rm *_out_p.dat