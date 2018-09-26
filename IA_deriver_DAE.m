%% Set up DAE
syms x1 x2 x3 x4 xM y1 y2 y3 y4 yM x1dot x2dot x3dot x4dot xMdot...
    y1dot y2dot y3dot y4dot yMdot x1ddot x2ddot x3ddot x4ddot xMddot y1ddot y2ddot y3ddot y4ddot yMddot real
syms th1 th2 th3 th4 th1dot th2dot th3dot th4dot th1ddot th2ddot th3ddot th4ddot real
syms N1x N1y N2x N2y N3x N3y N4x N4y N5x N5y real
syms m1 m2 m3 m4 M I1 I2 I3 I4 l1 l2 l3 l4 d1 d2 d3 d4 g real
syms a11 a12 a22 a23 a33 a34 real % muscle attachment points (a23 is muscle 2's attachment to bone 3)
syms M1 M2 M3 real % muscle force (scalars)

% arbitrary unit vectors, transformation matrices
i = [1 0 0]'; j = [0 1 0]'; k = [0 0 1]';

trans1 = [cos(th1) -sin(th1) 0; sin(th1) cos(th1) 0; 0 0 1];
trans2 = [cos(th2) -sin(th2) 0; sin(th2) cos(th2) 0; 0 0 1];
trans3 = [cos(th3) -sin(th3) 0; sin(th3) cos(th3) 0; 0 0 1];
trans4 = [cos(th4) -sin(th4) 0; sin(th4) cos(th4) 0; 0 0 1];

% vectors to each hinge
r_J1 = trans1*(l1*j); 
r_J2_J1 = trans2*(l2*j);
r_J2 = r_J1+r_J2_J1;
r_J3_J2 = trans3*(l3*j);
r_J3 = r_J2+r_J3_J2;
r_M_J3 = trans4*(l4*j);
r_M = r_M_J3+r_J3;

% vectors to each COM (from fixed frame)
r_g1 = trans1*(d1*j);
r_g2_J1 = trans2*(d2*j); r_g2 = r_J1 + r_g2_J1;
r_g3_J2 = trans3*(d3*j); r_g3 = r_J2 + r_g3_J2;
r_g4_J3 = trans4*(d4*j); r_g4 = r_J3 + r_g4_J3;

% muscle vectors
r_a11 = trans1*(a11*j);
r_a12_J1 = trans2*(a12*j); r_a12 = r_a12_J1+r_J1;

m11 = -(r_a11-r_a12)/(norm(r_a11-r_a12)); m12 = -m11; % m11 = muscle one on bone 1 direction, such that positive force value is a muscle contraction

r_a22_J1 = trans2*(a22*j); r_a22 = r_a22_J1+r_J1;
r_a23_J2 = trans3*(a23*j); r_a23 = r_a23_J2+r_J2;

m22 = -(r_a22-r_a23)/(norm(r_a22-r_a23)); m23 = -m22;

r_a33_J2 = trans3*(a33*j); r_a33 = r_a33_J2+r_J2;
r_a34_J3 = trans4*(a34*j); r_a34 = r_a34_J3+r_J3;

m33 = -(r_a33-r_a34)/(norm(r_a33-r_a34)); m34 = -m33;

%% LMB equations
% link 1
eq1 = N1x-N2x+dot(M1*m11,i)-m1*x1ddot;
eq2 = N1y-N2y+dot(M1*m11,j)+m1*g-m1*y1ddot;

% link 2
eq4 = N2x-N3x+dot(M1*m12,i)+dot(M2*m22,i)-m2*x2ddot;
eq5 = N2y-N3y+m2*g+dot(M1*m12,j)+dot(M2*m22,j)-m2*y2ddot;

% link 3
eq7 = N3x-N4x+dot(M2*m23,i)+dot(M3*m33,i)-m3*x3ddot;
eq8 = N3y-N4y+m3*g+dot(M2*m23,j)+dot(M3*m33,j)-m3*y3ddot;

% link 4
eq10 = N4x-N5x+dot(M3*m34,i)-m4*x4ddot;
eq11 = N4y-N5y+m4*g+dot(M3*m34,j)-m4*y4ddot;

% mass M
eq13 = N5x-M*xMddot;
eq14 = N5y+M*g-M*yMddot;

%% AMB
% link 1
eq3 = dot(I1*th1ddot*k-(cross(-d1*trans1*j,N1y*j)+cross(-d1*trans1*j,N1x*i)+cross((l1-d1)*trans1*j,-N2x*i)...
    +cross((l1-d1)*trans1*j,-N2y*j)+cross(r_a11-r_g1,M1*m11)),k);

% link 2
eq6 = dot(I2*th2ddot*k-(cross(-d2*trans2*j,N2x*i)+cross(-d2*trans2*j,N2y*j)+cross((l2-d2)*trans2*j,-N3x*i)...
    +cross((l2-d2)*trans2*j,-N3y*j)+cross(r_a12-r_g2,M1*m12)+cross(r_a22-r_g2,M2*m22)),k);

% link 3
eq9 = dot(I3*th3ddot*k-(cross(-d3*trans3*j,N3x*i)+cross(-d3*trans3*j,N3y*j)+cross((l3-d3)*trans3*j,-N4x*i)...
    +cross((l3-d3)*trans3*j,-N4y*j)+cross(r_a23-r_g3,M3*m23)+cross(r_a33-r_g3,M3*m33)),k);

% link 4
eq12 = dot(I4*th4ddot*k-(cross(-d4*trans4*j,N4x*i)+cross(-d4*trans4*j,N4y*j)+cross((l4-d4)*trans4*j,-N5x*i)...
    +cross((l4-d4)*trans4*j,-N5y*j)+cross(r_a34-r_g4,M3*m34)),k);

%% Constraint equations
% find acceleration of each COM
v_g1 = jacobian(r_g1,th1)*th1dot;
a_g1 = jacobian(v_g1,[th1 th1dot])*[th1dot th1ddot]';

v_g2 = jacobian(r_g2,[th1 th2])*[th1dot th2dot]';
a_g2 = jacobian(v_g2,[th1 th1dot th2 th2dot])*[th1dot th1ddot th2dot th2ddot]';

v_g3 = jacobian(r_g3,[th1 th2 th3])*[th1dot th2dot th3dot]';
a_g3 = jacobian(v_g3,[th1 th1dot th2 th2dot th3 th3dot])*[th1dot th1ddot th2dot th2ddot th3dot th3ddot]';

v_g4 = jacobian(r_g4,[th1 th2 th3 th4])*[th1dot th2dot th3dot th4dot]';
a_g4 = jacobian(v_g4,[th1 th1dot th2 th2dot th3 th3dot th4 th4dot])*[th1dot th1ddot th2dot th2ddot th3dot th3ddot th4dot th4ddot]';

v_M = jacobian(r_M,[th1 th2 th3 th4])*[th1dot th2dot th3dot th4dot]';
a_M = jacobian(v_M,[th1 th1dot th2 th2dot th3 th3dot th4 th4dot])*[th1dot th1ddot th2dot th2ddot th3dot th3ddot th4dot th4ddot]';

% constraints
eq15 = dot(a_g1-x1ddot*i-y1ddot*j,i);
eq16 = dot(a_g1-x1ddot*i-y1ddot*j,j);

eq17 = dot(a_g2-x2ddot*i-y2ddot*j,i);
eq18 = dot(a_g2-x2ddot*i-y2ddot*j,j);

eq19 = dot(a_g3-x3ddot*i-y3ddot*j,i);
eq20 = dot(a_g3-x3ddot*i-y3ddot*j,j);

eq21 = dot(a_g4-x4ddot*i-y4ddot*j,i);
eq22 = dot(a_g4-x4ddot*i-y4ddot*j,j);

eq23 = dot(a_M-xMddot*i-yMddot*j,i);
eq24 = dot(a_M-xMddot*i-yMddot*j,j);

%% Put equations in matrix, make function

[A,B] = equationsToMatrix([eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12 eq13 eq14 eq15 eq16 eq17 eq18 eq19 eq20 eq21 eq22 eq23 eq24],...
    [x1ddot x2ddot x3ddot x4ddot xMddot y1ddot y2ddot y3ddot y4ddot yMddot th1ddot th2ddot th3ddot th4ddot N1x N2x N3x N4x N5x N1y N2y N3y N4y N5y]);

matlabFunction(A,B,'file','IA_accel')

% x = A\B;