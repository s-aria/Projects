Point(1) = {0, 0.1, 0, 0.1};
Point(2) = {0, 0.4, 0, 0.1};
Point(3) = {1, 0.5, 0, 0.1};
Point(4) = {1, 0, 0, 0.1};
Point(5) = {0.5, 0, 0, 0.1};
Point(6) = {0.5, 0.1, 0, 0.1};
Line(1) = {5, 6};
Line(2) = {6, 1};
Line(3) = {1, 2};
Delete {
  Point{2};
}
Delete {
  Point{2};
}
Delete {
  Line{3, 2};
}
Delete {
  Point{2};
}
Point(7) = {0, 0.5, 0, 0.1};
Line(2) = {6, 1};
Line(3) = {1, 7};
Line(4) = {3, 3};
Line(5) = {4, 3};
Line(6) = {3, 7};
Line(7) = {5, 4};
Line Loop(8) = {6, -3, -2, -1, 7, 5};
Plane Surface(9) = {8};
