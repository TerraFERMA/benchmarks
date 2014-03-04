l=160;
p=1.0;
Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1}; 
}
Extrude {0, 0.5, 0} {
  Line{1};
}
Transfinite Line{1} = l+1;
Transfinite Line{2} = l+1;
Transfinite Line{-3} = (l+2)/2 Using Progression p;
Transfinite Line{-4} = (l+2)/2 Using Progression p;
Transfinite Surface{5} = {1,2,3,4} Right;
Extrude {0, 0.5, 0} {
  Line{2};
}
Transfinite Line{6} = l+1;
Transfinite Line{7} = (l+2)/2 Using Progression p;
Transfinite Line{8} = (l+2)/2 Using Progression p;
Transfinite Surface{9} = {3,4,5,6} Right;
// Bottom
Physical Line(3) = {1};
// Right
Physical Line(2) = {4, 8};
// Top
Physical Line(4) = {6};
// Left
Physical Line(1) = {7, 3};
// Internal
Physical Line(5) = {2};
Physical Surface(0) = {5, 9};
