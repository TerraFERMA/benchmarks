lz=${ncells};
p=${factor};
ax=${xaspect};
ay=${yaspect};
lx=lz*ax;
ly=lz*ay;
Point(1) = {0, 0, 0};
Extrude {ax, 0, 0} {
  Point{1}; 
}
Extrude {0, ay, 0} {
  Line{1};
}
Extrude {0, 0, 1.0} {
  Surface{5};
}
Transfinite Line{1} = lx+1 Using Bump p;
Transfinite Line{2} = lx+1 Using Bump p;
Transfinite Line{3} = ly+1 Using Bump p;
Transfinite Line{4} = ly+1 Using Bump p;
Transfinite Line{7} = lx+1 Using Bump p;
Transfinite Line{9} = lx+1 Using Bump p;
Transfinite Line{8} = ly+1 Using Bump p;
Transfinite Line{10} = ly+1 Using Bump p;
Transfinite Line{12} = lz+1 Using Bump p;
Transfinite Line{13} = lz+1 Using Bump p;
Transfinite Line{17} = lz+1 Using Bump p;
Transfinite Line{21} = lz+1 Using Bump p;
Transfinite Surface '*';
Transfinite Volume '*';
// Left
Physical Surface(1) = {14};
// Right
Physical Surface(2) = {22};
// Base
Physical Surface(3) = {5};
// Top
Physical Surface(4) = {27};
// Front
Physical Surface(5) = {18};
// Back
Physical Surface(6) = {26};
Physical Volume(0) = {1};
