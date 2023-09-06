# LoopSubdivision
To compile on linux 
module add qt/5.13.0
qmake -project QT+=opengl
qmake
make

./LoopSubdivisionRelease ../path_to/model.diredgenormal

Function:define a limit surface from a mesh by subdividing the polygons in the mesh.
