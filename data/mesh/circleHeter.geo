
// —————————CIRCLE———————————
// 	      3
//	  2   |   4
//     1 ——— 0 ——— 5 
//	  8   |   6
// 	      7
// input: (x0, y0) and r

Function CircleTemplate1
  // —
  x1 = x0-r;	y1 = y0;
  x5 = x0+r;	y5 = y0;
  x7 = x0;	y7 = y0-r;
  x3 = x0;	y3 = y0+r;
  // —
  dr = r/Sqrt(2);	
  x2 = x0-dr;	y2 = y0+dr;
  x8 = x0-dr;	y8 = y0-dr;
  x4 = x0+dr;	y4 = y0+dr;
  x6 = x0+dr;	y6 = y0-dr;
  // points
  p0 = newp; Point(p0) = { x0, y0, 0, cl }; circlePoint1[9*cc+0] = p0;
  p1 = newp; Point(p1) = { x1, y1, 0, cl }; circlePoint1[9*cc+1] = p1;
  p2 = newp; Point(p2) = { x2, y2, 0, cl }; circlePoint1[9*cc+2] = p2;
  p3 = newp; Point(p3) = { x3, y3, 0, cl }; circlePoint1[9*cc+3] = p3;
  p4 = newp; Point(p4) = { x4, y4, 0, cl }; circlePoint1[9*cc+4] = p4;
  p5 = newp; Point(p5) = { x5, y5, 0, cl }; circlePoint1[9*cc+5] = p5;
  p6 = newp; Point(p6) = { x6, y6, 0, cl }; circlePoint1[9*cc+6] = p6;
  p7 = newp; Point(p7) = { x7, y7, 0, cl }; circlePoint1[9*cc+7] = p7;
  p8 = newp; Point(p8) = { x8, y8, 0, cl }; circlePoint1[9*cc+8] = p8;
  // lines
  l1 = newl; Circle(l1) = { p1, p0, p2 }; circleLines1[8*cc+0] = l1;
  l2 = newl; Circle(l2) = { p2, p0, p3 }; circleLines1[8*cc+1] = l2;
  l3 = newl; Circle(l3) = { p3, p0, p4 }; circleLines1[8*cc+2] = l3;
  l4 = newl; Circle(l4) = { p4, p0, p5 }; circleLines1[8*cc+3] = l4;
  l5 = newl; Circle(l5) = { p5, p0, p6 }; circleLines1[8*cc+4] = l5;
  l6 = newl; Circle(l6) = { p6, p0, p7 }; circleLines1[8*cc+5] = l6;
  l7 = newl; Circle(l7) = { p7, p0, p8 }; circleLines1[8*cc+6] = l7;
  l8 = newl; Circle(l8) = { p8, p0, p1 }; circleLines1[8*cc+7] = l8;
  l11 = newl; Line(l11) = { p1, p0};
  l12 = newl; Line(l12) = { p0, p5};
  l13 = newl; Line(l13) = { p3, p0};
  l14 = newl; Line(l14) = { p0, p7};
  //  line loop
  ll1 = newll; Line Loop(ll1) = { l1, l2,  l13, -l11 };
  ll2 = newll; Line Loop(ll2) = { l3, l4, -l12, -l13 };
  ll3 = newll; Line Loop(ll3) = { l5, l6, -l14,  l12 };
  ll4 = newll; Line Loop(ll4) = { l7, l8,  l11,  l14 };
  s1 = news; Ruled Surface(s1) = { ll1 };
  s2 = news; Ruled Surface(s2) = { ll2 };
  s3 = news; Ruled Surface(s3) = { ll3 };
  s4 = news; Ruled Surface(s4) = { ll4 };
  llc = newll; Line Loop(llc) = { l1, l2, l3, l4, l5, l6, l7, l8 };
  circleLineLoop1[cc] = llc;
  circleSurf1[4*cc+0] = s1;
  circleSurf1[4*cc+1] = s2;
  circleSurf1[4*cc+2] = s3;
  circleSurf1[4*cc+3] = s4;

  l101 = newl; Line(l101) = {p0, p2}; Line{l101} In Surface{s1}; 
  l102 = newl; Line(l102) = {p0, p4}; Line{l102} In Surface{s2}; 
  l103 = newl; Line(l103) = {p0, p6}; Line{l103} In Surface{s3}; 
  l104 = newl; Line(l104) = {p0, p8}; Line{l104} In Surface{s4}; 
Return


Function CircleTemplate2
  // —
  x1 = x0-r;	y1 = y0;
  x5 = x0+r;	y5 = y0;
  x7 = x0;	y7 = y0-r;
  x3 = x0;	y3 = y0+r;
  // —
  dr = r/Sqrt(2);	
  x2 = x0-dr;	y2 = y0+dr;
  x8 = x0-dr;	y8 = y0-dr;
  x4 = x0+dr;	y4 = y0+dr;
  x6 = x0+dr;	y6 = y0-dr;
  // points
  p0 = newp; Point(p0) = { x0, y0, 0, cl }; circlePoint2[9*cc+0] = p0;
  p1 = newp; Point(p1) = { x1, y1, 0, cl }; circlePoint2[9*cc+1] = p1;
  p2 = newp; Point(p2) = { x2, y2, 0, cl }; circlePoint2[9*cc+2] = p2;
  p3 = newp; Point(p3) = { x3, y3, 0, cl }; circlePoint2[9*cc+3] = p3;
  p4 = newp; Point(p4) = { x4, y4, 0, cl }; circlePoint2[9*cc+4] = p4;
  p5 = newp; Point(p5) = { x5, y5, 0, cl }; circlePoint2[9*cc+5] = p5;
  p6 = newp; Point(p6) = { x6, y6, 0, cl }; circlePoint2[9*cc+6] = p6;
  p7 = newp; Point(p7) = { x7, y7, 0, cl }; circlePoint2[9*cc+7] = p7;
  p8 = newp; Point(p8) = { x8, y8, 0, cl }; circlePoint2[9*cc+8] = p8;
  // lines
  l1 = newl; Circle(l1) = { p1, p0, p2 }; circleLines2[8*cc+0] = l1;
  l2 = newl; Circle(l2) = { p2, p0, p3 }; circleLines2[8*cc+1] = l2;
  l3 = newl; Circle(l3) = { p3, p0, p4 }; circleLines2[8*cc+2] = l3;
  l4 = newl; Circle(l4) = { p4, p0, p5 }; circleLines2[8*cc+3] = l4;
  l5 = newl; Circle(l5) = { p5, p0, p6 }; circleLines2[8*cc+4] = l5;
  l6 = newl; Circle(l6) = { p6, p0, p7 }; circleLines2[8*cc+5] = l6;
  l7 = newl; Circle(l7) = { p7, p0, p8 }; circleLines2[8*cc+6] = l7;
  l8 = newl; Circle(l8) = { p8, p0, p1 }; circleLines2[8*cc+7] = l8;
  l11 = newl; Line(l11) = { p1, p0};
  l12 = newl; Line(l12) = { p0, p5};
  l13 = newl; Line(l13) = { p3, p0};
  l14 = newl; Line(l14) = { p0, p7};
  //  line loop
  ll1 = newll; Line Loop(ll1) = { l1, l2,  l13, -l11 };
  ll2 = newll; Line Loop(ll2) = { l3, l4, -l12, -l13 };
  ll3 = newll; Line Loop(ll3) = { l5, l6, -l14,  l12 };
  ll4 = newll; Line Loop(ll4) = { l7, l8,  l11,  l14 };
  s1 = news; Ruled Surface(s1) = { ll1 };
  s2 = news; Ruled Surface(s2) = { ll2 };
  s3 = news; Ruled Surface(s3) = { ll3 };
  s4 = news; Ruled Surface(s4) = { ll4 };
  llc = newll; Line Loop(llc) = { l1, l2, l3, l4, l5, l6, l7, l8 };
  circleLineLoop2[cc] = llc;
  circleSurf2[4*cc+0] = s1;
  circleSurf2[4*cc+1] = s2;
  circleSurf2[4*cc+2] = s3;
  circleSurf2[4*cc+3] = s4;

  l101 = newl; Line(l101) = {p0, p2}; Line{l101} In Surface{s1}; 
  l102 = newl; Line(l102) = {p0, p4}; Line{l102} In Surface{s2}; 
  l103 = newl; Line(l103) = {p0, p6}; Line{l103} In Surface{s3}; 
  l104 = newl; Line(l104) = {p0, p8}; Line{l104} In Surface{s4}; 
Return



Function CircleTemplate3
  // —
  x1 = x0-r;	y1 = y0;
  x5 = x0+r;	y5 = y0;
  x7 = x0;	y7 = y0-r;
  x3 = x0;	y3 = y0+r;
  // —
  dr = r/Sqrt(2);	
  x2 = x0-dr;	y2 = y0+dr;
  x8 = x0-dr;	y8 = y0-dr;
  x4 = x0+dr;	y4 = y0+dr;
  x6 = x0+dr;	y6 = y0-dr;
  // points
  p0 = newp; Point(p0) = { x0, y0, 0, cl }; circlePoint3[9*cc+0] = p0;
  p1 = newp; Point(p1) = { x1, y1, 0, cl }; circlePoint3[9*cc+1] = p1;
  p2 = newp; Point(p2) = { x2, y2, 0, cl }; circlePoint3[9*cc+2] = p2;
  p3 = newp; Point(p3) = { x3, y3, 0, cl }; circlePoint3[9*cc+3] = p3;
  p4 = newp; Point(p4) = { x4, y4, 0, cl }; circlePoint3[9*cc+4] = p4;
  p5 = newp; Point(p5) = { x5, y5, 0, cl }; circlePoint3[9*cc+5] = p5;
  p6 = newp; Point(p6) = { x6, y6, 0, cl }; circlePoint3[9*cc+6] = p6;
  p7 = newp; Point(p7) = { x7, y7, 0, cl }; circlePoint3[9*cc+7] = p7;
  p8 = newp; Point(p8) = { x8, y8, 0, cl }; circlePoint3[9*cc+8] = p8;
  // lines
  l1 = newl; Circle(l1) = { p1, p0, p2 }; circleLines3[8*cc+0] = l1;
  l2 = newl; Circle(l2) = { p2, p0, p3 }; circleLines3[8*cc+1] = l2;
  l3 = newl; Circle(l3) = { p3, p0, p4 }; circleLines3[8*cc+2] = l3;
  l4 = newl; Circle(l4) = { p4, p0, p5 }; circleLines3[8*cc+3] = l4;
  l5 = newl; Circle(l5) = { p5, p0, p6 }; circleLines3[8*cc+4] = l5;
  l6 = newl; Circle(l6) = { p6, p0, p7 }; circleLines3[8*cc+5] = l6;
  l7 = newl; Circle(l7) = { p7, p0, p8 }; circleLines3[8*cc+6] = l7;
  l8 = newl; Circle(l8) = { p8, p0, p1 }; circleLines3[8*cc+7] = l8;
  l11 = newl; Line(l11) = { p1, p0};
  l12 = newl; Line(l12) = { p0, p5};
  l13 = newl; Line(l13) = { p3, p0};
  l14 = newl; Line(l14) = { p0, p7};
  //  line loop
  ll1 = newll; Line Loop(ll1) = { l1, l2,  l13, -l11 };
  ll2 = newll; Line Loop(ll2) = { l3, l4, -l12, -l13 };
  ll3 = newll; Line Loop(ll3) = { l5, l6, -l14,  l12 };
  ll4 = newll; Line Loop(ll4) = { l7, l8,  l11,  l14 };
  s1 = news; Ruled Surface(s1) = { ll1 };
  s2 = news; Ruled Surface(s2) = { ll2 };
  s3 = news; Ruled Surface(s3) = { ll3 };
  s4 = news; Ruled Surface(s4) = { ll4 };
  llc = newll; Line Loop(llc) = { l1, l2, l3, l4, l5, l6, l7, l8 };
  circleLineLoop3[cc] = llc;
  circleSurf3[4*cc+0] = s1;
  circleSurf3[4*cc+1] = s2;
  circleSurf3[4*cc+2] = s3;
  circleSurf3[4*cc+3] = s4;

  l101 = newl; Line(l101) = {p0, p2}; Line{l101} In Surface{s1}; 
  l102 = newl; Line(l102) = {p0, p4}; Line{l102} In Surface{s2}; 
  l103 = newl; Line(l103) = {p0, p6}; Line{l103} In Surface{s3}; 
  l104 = newl; Line(l104) = {p0, p8}; Line{l104} In Surface{s4}; 
Return


