Point(0) = {0,0,-1};
Point(1) = {1,0,-1};
Point(2) = {0,1,-1};

Point(3) = {0,0,1};
Point(4) = {1,0,1};
Point(5) = {0,1,1};

Line(6) = {0,1};
Line(7) = {1,2};
Line(8) = {2,0};

Line(9) = {0,3};
Line(10) = {1,4};
Line(11) = {2,5};

Line(12) = {3,4};
Line(13) = {4,5};
Line(14) = {5,3};


Line Loop(15) = {6,7,8};
Line Loop(16) = {6,10,-12,-9};
Line Loop(17) = {7,11,-13,-10};
Line Loop(18) = {8,9,-14,-11};
Line Loop(19) = {12,13,14};

Surface(15) = {15};
Surface(16) = {16};
Surface(17) = {17};
Surface(18) = {18};
Surface(19) = {19};


Surface Loop(20) = {15,16,17,18,19};
Volume(20) = {20};


Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Mesh 3;