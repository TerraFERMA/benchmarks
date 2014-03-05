inner = 0.000775;
third = inner;
twothird = 0.0016;
top = 0.0024;
middle = 0.005;
outer = 0.012;
Point(1) = {0.009, 0, 0, inner};
Point(2) = {0.009, 0.003, 0, inner};
Point(5) = {0, 0.003, 0, inner};
Point(6) = {0.085, 0, 0, outer};
Point(19) = {0.085, 0.05, 0, outer};
Point(20) = {0.085, 0.1, 0, outer};
Point(9) = {0.085, 0.161, 0, outer};
Point(10) = {0, 0.161, 0, top};
Point(11) = {0, 0.05, 0, third};
Point(12) = {0, 0.1, 0, twothird};
Point(13) = {0.02, 0, 0, inner};
Point(14) = {0.02, 0.05, 0, third};
Point(15) = {0.02, 0.1, 0, twothird};
Point(22) = {0.02, 0.161, 0, top};
Point(16) = {0.05, 0, 0, middle};
Point(17) = {0.05, 0.05, 0, middle};
Point(18) = {0.05, 0.1, 0, middle};
Point(21) = {0.05, 0.161, 0, middle};
Line(1) = {5, 11};
Line(2) = {11, 12};
Line(3) = {12, 10};
Line(4) = {10, 22};
Line(5) = {22, 21};
Line(6) = {21, 9};
Line(7) = {9, 20};
Line(8) = {20, 19};
Line(9) = {19, 6};
Line(10) = {6, 16};
Line(11) = {16, 13};
Line(12) = {13, 1};
Line(13) = {1, 2};
Line(14) = {2, 5};
Line(15) = {13, 14};
Line(16) = {14, 15};
Line(17) = {15, 22};
Line(18) = {21, 18};
Line(19) = {18, 17};
Line(20) = {17, 16};
Line Loop(21) = {2, 3, 4, -17, -16, -15, 12, 13, 14, 1};
Plane Surface(22) = {21};
Line Loop(23) = {16, 17, 5, 18, 19, 20, 11, 15};
Plane Surface(24) = {23};
Line Loop(25) = {20, -10, -9, -8, -7, -6, 18, 19};
Plane Surface(26) = {25};
Physical Line(27) = {2, 3, 1};
Physical Line(28) = {4, 5, 6};
Physical Line(29) = {7, 8, 9};
Physical Line(30) = {10, 11, 12};
Physical Surface(31) = {22};
Physical Surface(32) = {24};
Physical Surface(33) = {26};
Physical Line(34) = {13, 14};
