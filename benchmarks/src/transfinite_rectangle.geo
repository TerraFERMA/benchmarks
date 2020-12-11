ly=${ncells};
pv=${factor};
ph=1.0;
a=${aspect};
lx=ly*a;
Point(1) = {0, 0, 0};
Extrude {a, 0, 0} {
  Point{1}; 
}
Extrude {0, 1.0, 0} {
  Line{1};
}
Transfinite Line{1} = lx+1 Using Bump ph;
Transfinite Line{2} = lx+1 Using Bump ph;
Transfinite Line{3} = ly+1 Using Bump pv;
Transfinite Line{4} = ly+1 Using Bump pv;
Transfinite Surface{5} = {1,2,3,4} Right;
// Bottom
Physical Line(3) = {1};
// Right
Physical Line(2) = {4};
// Top
Physical Line(4) = {2};
// Left
Physical Line(1) = {3};
Physical Surface(0) = {5};
