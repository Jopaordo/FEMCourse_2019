Point(0) = {-1,-1,-1};
Point(1) = {1,-1,-1};
Point(2) = {1,1,-1};
Point(3) = {-1,1,-1};

Point(4) = {-1,-1,1};
Point(5) = {1,-1,1};
Point(6) = {1,1,1};
Point(7) = {-1,1,1};


Line(8) = {0,1};
Line(9) = {1,2};
Line(10) = {2,3};
Line(11) = {3,0};

Line(12) = {0,4};
Line(13) = {1,5};
Line(14) = {2,6};
Line(15) = {3,7};

Line(16) = {4,5};
Line(17) = {5,6};
Line(18) = {6,7};
Line(19) = {7,4};

Line Loop(20) = {8,9,10,11};
Line Loop(21) = {8,13,-16,-12};
Line Loop(22) = {9,14,-17,-13};
Line Loop(23) = {10,15,-18,-14};
Line Loop(24) = {11,12,-19,-15};
Line Loop(25) = {16,17,18,19};

Surface(20) = {20};
Surface(21) = {21};
Surface(22) = {22};
Surface(23) = {23};
Surface(24) = {24};
Surface(25) = {25};

Surface Loop(26) = {20,21,22,23,24,25};
Volume(26) = {26};


Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Mesh 3;