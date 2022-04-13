// Gmsh project created on Thu Apr 07 13:27:39 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {2, 1, 0, 1.0};
//+
Point(5) = {2, 0, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 1};
//+
Curve Loop(1) = {2, 3, 4, 1};

//+
Circle(5) = {0.5, 0.5, 0, 0.25, 0, 2*Pi};
//+
Curve Loop(2) = {2, 3, 4, 1};
//+
Curve Loop(3) = {5};
//+
Plane Surface(1) = {2, 3};
//+
Physical Curve("Inlet", 6) = {1};
//+
Physical Curve("Wall", 7) = {2, 4};
//+
Physical Curve("Outlet", 8) = {3};
//+
Physical Curve("Airfoil", 9) = {5};
//+
Physical Surface("Fluid", 10) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {2, 4} = 30 Using Progression 1;
//+
Transfinite Curve {1, 3} = 15 Using Progression 1;
//+
Recombine Surface {1};
