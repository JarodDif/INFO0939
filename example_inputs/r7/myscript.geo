Mesh.SurfaceEdges = 0;
Mesh.SurfaceFaces = 0;
Mesh.VolumeEdges = 0;
For i In {0:PostProcessing.NbViews-1}
  View[i].RangeType = 2; // custom range
  View[i].CustomMin = -0.5;
  View[i].CustomMax = 0.5;
  View[i].RaiseZ = 3;
  If (i == 0)
    View[i].Name = "R7";
  Else
    View[i].ShowScale = 0;
  EndIf
EndFor
