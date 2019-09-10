Point(0) = {-1,-1,0};
Point(1) = {1,-1,0};
Point(2) = {1,1,0};
Point(3) = {-1,1,0};
Point(4) = {0,0,1};

Line(5) = {0,1};
Line(6) = {1,2};
Line(7) = {2,3};
Line(8) = {3,0};

Line(9) = {0,4};
Line(10) = {1,4};
Line(11) = {2,4};
Line(12) = {3,4};

Line Loop(13) = {5,6,7,8};
Line Loop(14) = {5,10,-9};
Line Loop(15) = {6,11,-10};
Line Loop(16) = {7,12,-11};
Line Loop(17) = {8,9,-12};

Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};
Plane Surface(17) = {17};

Surface Loop(18) = {13,14,15,16,17};
Volume(18) = {18};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

// Mesh 3;