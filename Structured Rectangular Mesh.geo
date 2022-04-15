
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {2, 1, 0, 1.0};
//+
Point(4) = {2, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Physical Curve("Wall", 5) = {2, 4};
//+
Physical Curve("Inlet", 6) = {1};
//+
Physical Curve("Outlet", 7) = {3};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 3} = 5 Using Progression 1;
//+
Transfinite Curve {2, 4} = 10 Using Progression 1;
//+
Recombine Surface {1};
//+
Physical Surface("Fluid", 8) = {1};
