function [A,B] = IA_accel(I1,I2,I3,I4,M,M1,M2,M3,a11,a12,a22,a23,a33,a34,d1,d2,d3,d4,g,l1,l2,l3,l4,m1,m2,m3,m4,th1,th2,th3,th4,th1dot,th2dot,th3dot,th4dot)
%IA_ACCEL
%    [A,B] = IA_ACCEL(I1,I2,I3,I4,M,M1,M2,M3,A11,A12,A22,A23,A33,A34,D1,D2,D3,D4,G,L1,L2,L3,L4,M1,M2,M3,M4,TH1,TH2,TH3,TH4,TH1DOT,TH2DOT,TH3DOT,TH4DOT)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    23-Sep-2018 16:08:12

t2 = cos(th1);
t3 = sin(th1);
t4 = d1-l1;
t5 = cos(th2);
t6 = sin(th2);
t7 = d2-l2;
t8 = cos(th3);
t9 = sin(th3);
t10 = d3-l3;
t11 = cos(th4);
t12 = sin(th4);
t13 = d4-l4;
A = reshape([-m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-M,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,-m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-M,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,I1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-d1.*t2,-d1.*t3,-l1.*t2,-l1.*t3,-l1.*t2,-l1.*t3,-l1.*t2,-l1.*t3,-l1.*t2,-l1.*t3,0.0,0.0,0.0,0.0,0.0,I2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-d2.*t5,-d2.*t6,-l2.*t5,-l2.*t6,-l2.*t5,-l2.*t6,-l2.*t5,-l2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,I3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-d3.*t8,-d3.*t9,-l3.*t8,-l3.*t9,-l3.*t8,-l3.*t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,I4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-d4.*t11,-d4.*t12,-l4.*t11,-l4.*t12,1.0,0.0,-d1.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t2.*t4,1.0,0.0,-d2.*t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t5.*t7,1.0,0.0,-d3.*t8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t8.*t10,1.0,0.0,-d4.*t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t11.*t13,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,-d1.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t3.*t4,0.0,1.0,-d2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t6.*t7,0.0,1.0,-d3.*t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t9.*t10,0.0,1.0,-d4.*t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t12.*t13,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[24,24]);
if nargout > 1
    t16 = a11.*t3;
    t17 = a12.*t6;
    t18 = l1.*t3;
    t19 = -t16+t17+t18;
    t14 = abs(t19);
    t21 = a11.*t2;
    t22 = a12.*t5;
    t23 = l1.*t2;
    t24 = -t21+t22+t23;
    t15 = abs(t24);
    t20 = t14.^2;
    t25 = t15.^2;
    t26 = t20+t25;
    t27 = 1.0./sqrt(t26);
    t28 = M1.*t19.*t27;
    t31 = a22.*t6;
    t32 = a23.*t9;
    t33 = l2.*t6;
    t34 = -t31+t32+t33;
    t29 = abs(t34);
    t36 = a22.*t5;
    t37 = a23.*t8;
    t38 = l2.*t5;
    t39 = -t36+t37+t38;
    t30 = abs(t39);
    t35 = t29.^2;
    t40 = t30.^2;
    t41 = t35+t40;
    t42 = 1.0./sqrt(t41);
    t43 = M2.*t34.*t42;
    t46 = a33.*t9;
    t47 = a34.*t12;
    t48 = l3.*t9;
    t49 = -t46+t47+t48;
    t44 = abs(t49);
    t51 = a33.*t8;
    t52 = a34.*t11;
    t53 = l3.*t8;
    t54 = -t51+t52+t53;
    t45 = abs(t54);
    t50 = t44.^2;
    t55 = t45.^2;
    t56 = t50+t55;
    t57 = 1.0./sqrt(t56);
    t58 = M3.*t49.*t57;
    t59 = th1dot.^2;
    t60 = th2dot.^2;
    t61 = th3dot.^2;
    t62 = l1.*t2.*t59;
    t63 = th4dot.^2;
    t64 = l2.*t5.*t60;
    t65 = l3.*t8.*t61;
    B = [t28;-g.*m1-M1.*t24.*t27;M1.*t19.*t27.*(t21-d1.*t2)-M1.*t24.*t27.*(t16-d1.*t3);-t28+t43;-g.*m2+M1.*t24.*t27-M2.*t39.*t42;-M1.*t19.*t27.*(t22-d2.*t5)+M1.*t24.*t27.*(t17-d2.*t6)+M2.*t34.*t42.*(t36-d2.*t5)-M2.*t39.*t42.*(t31-d2.*t6);-t43+t58;-g.*m3+M2.*t39.*t42-M3.*t54.*t57;-M3.*t34.*t42.*(t37-d3.*t8)+M3.*t39.*t42.*(t32-d3.*t9)+M3.*t49.*t57.*(t51-d3.*t8)-M3.*t54.*t57.*(t46-d3.*t9);-t58;-g.*m4+M3.*t54.*t57;-M3.*t49.*t57.*(t52-d4.*t11)+M3.*t54.*t57.*(t47-d4.*t12);0.0;-M.*g;-d1.*t3.*t59;d1.*t2.*t59;-d2.*t6.*t60-l1.*t3.*t59;t62+d2.*t5.*t60;-d3.*t9.*t61-l1.*t3.*t59-l2.*t6.*t60;t62+t64+d3.*t8.*t61;-d4.*t12.*t63-l1.*t3.*t59-l2.*t6.*t60-l3.*t9.*t61;t62+t64+t65+d4.*t11.*t63;-l1.*t3.*t59-l2.*t6.*t60-l3.*t9.*t61-l4.*t12.*t63;t62+t64+t65+l4.*t11.*t63];
end
