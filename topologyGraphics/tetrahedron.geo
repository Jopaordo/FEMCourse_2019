Point(0) = {0,0,0};
Point(1) = {1,0,0};
Point(2) = {0,1,0};
Point(3) = {0,0,1};

Line(4) = {0,1};
Line(5) = {1,2};
Line(6) = {2,0};
Line(7) = {0,3};
Line(8) = {1,3};
Line(9) = {2,3};

Line Loop(10) = {4,5,6};
Line Loop(11) = {4,8,-7};
Line Loop(12) = {5,9,-8};
Line Loop(13) = {6,7,-9};

Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

Surface Loop(14) = {10,11,12,13};
Volume(14) = {14};

Physical Surface(2) = {10};
Physical Surface(3) = {11};