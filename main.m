read_obj
verts = ans;
DT = delaunayTriangulation(verts.vertices);
tetramesh(DT);