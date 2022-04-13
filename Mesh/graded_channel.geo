// Cylinder with radius 0.1 in 2.2x3.2 rectangular domain
// Cylinder centre at (0, 0)

SetFactory("OpenCASCADE");

lc1 = 1e-2;

Point(1) = {   0,   0, 0 }; // bottom left corner
Point(2) = {   0, 0.2, 0 }; // top left corner
Point(3) = { 2.0, 0.2, 0 }; // top right corner
Point(4) = { 2.0,   0, 0 }; // bottom right corner

// external boundaries
Line(1)   = {1,2}; // inlet
Line(2)   = {2,3}; // top boundary
Line(3)   = {3,4}; // outlet
Line(4)   = {4,1}; // bottom boundary

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Line("inlet", 1)       = {1}; // inlet BC
Physical Line("outlet", 2)      = {3}; // outlet BC
Physical Line("bottom", 3)      = {4}; // bottom BC
Physical Line("top", 4)         = {2}; // top BC

Transfinite Line{1} = 20;
Transfinite Line{2} = 20;
Transfinite Line{3} = 100;
Transfinite Line{4} = 100;

Mesh.Algorithm = 8;

Transfinite Surface{1} = {1,2,3,4};

Recombine Surface {1};

Physical Surface(1) = {1};

