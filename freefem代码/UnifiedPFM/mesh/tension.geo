//设置内核
SetFactory("OpenCASCADE");
//定义点附近网格尺寸
lc = 0.5;
//定义切口尺寸
length = 0.5;
width = 0.01;
//定义点
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0.5+width/2, 0, lc};
Point(6) = {length, 0.5+width/2, 0, lc};
Point(7) = {length, 0.5-width/2, 0, lc};
Point(8) = {0, 0.5-width/2, 0, lc};
//定义曲线
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
//定义封闭环
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
//定义面
Plane Surface(1) = {1};
//划分物理区域
Physical Curve("bottom") = {1};
Physical Curve("up") = {3};
Physical Surface("plate") = {1};

//定义网格加密场
Field[1] = Box;
Field[1].Thickness = 0.5;  
Field[1].VIn = 0.005;
Field[1].VOut = 0.1;
Field[1].XMax = 1.1;
Field[1].XMin = 0.49;
Field[1].YMax = 0.52;
Field[1].YMin = 0.48;
Field[1].ZMax = 0.1;
Field[1].ZMin = -0.1;

// 应用网格场
Background Field = 1;
//绘制网格
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Smoothing = 5;
Mesh.Algorithm = 8;

Mesh 2;
//保存网格文件
Mesh.MshFileVersion = 2.2; // 设置为版本 2 ASCII
Mesh.Format = 1; // 设置为 ASCII 格式（1 = ASCII, 2 = Binary）
Save "tension.msh";
