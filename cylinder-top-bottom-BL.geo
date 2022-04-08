// Cylinder with radius 0.1 in 2.2x3.2 rectangular domain
// Cylinder centre at (0, 0)

SetFactory("OpenCASCADE");

lc1 = 1e-2;

Point(1) = {   0,   0 , 0, lc1}; // cylinder centre
Point(2) = {   0, -0.1, 0, lc1}; // cylinder bottom point
Point(3) = {   0,  0.1, 0, lc1}; // cylinder top point
Point(4) = {-1.1, -1.1, 0, 5*lc1};
Point(5) = {-1.1,  1.1, 0, 5*lc1};
Point(6) = { 2.1,  1.1, 0, 5*lc1};
Point(7) = { 2.1, -1.1, 0, 5*lc1};

// cylinder LE and TE
Circle(1) = {2,1,3}; // cylinder LE
Circle(2) = {3,1,2}; // cylinder TE

// external boundaries
Line(3)   = {4,5}; // inlet
Line(4)   = {5,6}; // top boundary
Line(5)   = {6,7}; // outlet
Line(6)   = {7,4}; // bottom boundary

Line Loop(1) = {3,4,5,6};
Line Loop(2) = {1,2};

Plane Surface(1) = {1,2};

Physical Line("inlet", 1)       = {3}; // inlet BC
Physical Line("outlet", 2)      = {5}; // outlet BC
Physical Line("bottom", 3)      = {6}; // bottom BC
Physical Line("top", 4)         = {4}; // top BC
Physical Line("cylinder", 5)    = {1,2}; // cylinder BC

Physical Surface("flow") = {1};

// cylinder boundary layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {1,2,4,6};
Field[1].hfar = 0.05;
Field[1].hwall_n = 0.005;
Field[1].thickness = 0.05;
Field[1].ratio = 1.1;
Field[1].IntersectMetrics = 1;
Field[1].Quads = 1;
BoundaryLayer Field = 1;

Transfinite Line{3,5} = 60 Using Bump 0.05;

// top boundary layer
// Field[2] = BoundaryLayer;
// Field[2].EdgesList = {4,6};
// Field[2].hfar = 0.05;
// Field[2].hwall_n = 0.005;
// Field[2].thickness = 0.05;
// Field[2].ratio = 1.1;
// Field[2].IntersectMetrics = 1;
// Field[2].Quads = 1;
// BoundaryLayer Field = 2;

// Extrusion from top boundary
// Extrude{0,0.05,0} {
//     Line{4}; Layers{5}; Recombine;}

// Extrusion from bottom boundary
// Extrude{0,-0.05,0} {
//     Line{6}; Layers{5}; Recombine;}

Mesh.Algorithm = 6;
Recombine Surface {1};


