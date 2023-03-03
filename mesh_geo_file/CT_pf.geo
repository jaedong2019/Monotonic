// Gmsh project created on Mon Apr 12 20:20:44 2020
SetFactory("OpenCASCADE");

m = 50; // dimensionless

f0 = 0.133;

a0 = 24.5/m;
a1 = 25/m;
h = 24./m;
L = 50./m;
eta = 0.001/m;
R = 5.2/m;
x = 10./m;
y = 11./m;

beta = 50*Pi/180;  //angle of contact surface
a1x = R-R*Cos(beta);
a1y = R*Sin(beta);

 
//+
Point(1) = {0, eta, 0, f0};
//+
Point(100) ={0, -eta, 0, f0};
//+
Point(2) = {0, h, 0, f0};
//+
Point(3) = {0, -h, 0, f0};
//+
Point(4) = {L, h, 0, f0};
//+
Point(5) = {L, -h, 0, f0};
//+
Point(6) = {a0, eta, 0, f0};
//+
Point(7) = {a0, -eta, 0, f0};
//+
Point(8) = {a1, 0, 0, f0};

//+
Point(9) = {x, y, 0, f0};
//+
Point(10) = {x-R+a1x, y+a1y, 0, f0};
//+
Point(11) = {x+R-a1x, y+a1y, 0, f0};
//+
Point(12) = {x, y-R, 0, f0};
//+
Point(13) = {x, y+R, 0, f0};
//+
Point(14) = {x, -y, 0, f0};
//+
Point(15) = {x, -y-R, 0, f0};
//+
Point(16) = {x, -y+R, 0, f0};
//+
Point(17) = {x-R+a1x, -y-a1y, 0, f0};
//+
Point(18) = {x+R-a1x, -y-a1y, 0, f0};




//+
Circle(1) = {13, 9, 11};
//+
Circle(2) = {11, 9, 12};
//+
Circle(3) = {12, 9, 10};
//+
Circle(4) = {10, 9, 13};
//+
Circle(5) = {16, 14, 18};
//+
Circle(6) = {18, 14, 15};
//+
Circle(7) = {15, 14, 17};
//+
Circle(8) = {17, 14, 16};
//+
Line(9) = {1, 6};
//+
Line(10) = {6, 8};
//+
Line(11) = {8, 7};
//+
Line(12) = {7, 100};
//+
Line(13) = {100, 3};
//+
Line(14) = {3, 5};
//+
Line(15) = {5, 4};
//+
Line(16) = {4, 2};
//+
Line(17) = {2, 1};
//+
Curve Loop(1) = {16, 17, 9, 10, 11, 12, 13, 14, 15};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Curve Loop(3) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2, 3};


//+
Physical Surface("1") = {1};
//+
Physical Curve("2") = {1, 4, 3, 2};
//+
Physical Curve("3") = {5, 6, 7, 8};
//+
Physical Curve("4") = {4, 1};

