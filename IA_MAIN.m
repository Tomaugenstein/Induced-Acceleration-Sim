function IA_MAIN

% Tom Augenstein, NeuRRo Lab

% The purpose of this code is to simulate and animate a test case of induced acceleration.
% The script references of MATLAB-generated function, "IA_accel.m",
% built from a symbolic function IA_deriver_DAE.m. This script uses
% Newton-Euler mechanics in maximal coordinates, also known as DAE, to
% solve for the 4 degrees of freedom of the system. The objective of this
% set of functions is to investigate whether or not all muscles in a system
% can create acceleration of the center of mass when activated
% individually.

hold on

% initialize constants, options
p.m1 = 1; p.m2 = 1; p.m3 = 1; p.m4 = 1; p.M = 20; % mass of each bone and top mass
mTot = p.m1+p.m2+p.m3+p.m4+p.M;

p.l1 = 1; p.l2 = 2; p.l3 = 2; p.l4 = 1; % bone lengths
p.d1 = p.l1/2; p.d2 = p.l2/2; p.d3 = p.l3/2; p.d4 = p.l4/2; % positions of center of mass on each bone
p.g = 0; % 9.81 = gravity on, 0 = gravity off
p.I1 = 1; p.I2 = 1; p.I3 = 1; p.I4 = 1; % bone moment of inertia
p.M1 = 22; p.M2 = 21; p.M3 = 20.00001; % muscles forces
p.a11 = 0; p.a12 = p.d2; p.a22 = p.d2; p.a23 = p.d3; p.a33 = p.d3; p.a34 = p.l4; % muscle attachment points (a23 = attachment point of muscle 2 on bone 3)

% intergration tolerances
opts.RelTol = 1e-6; opts.AbsTol = 1e-6;

% timeframe of integration
tf = 1;

% initial conditions of generalized coordinates, coordinates are angles of
% joints, measured from fixed frame (not from parent link)
Z_o = [pi/4 0 -pi/4 0 pi/4 0 -pi/4 0]'; % th1 th1dot th2 th2dot etc.

% numerically integrate (ode45 = runge-kutta 4 and runge-kutta 5)
[t, Z] = ode45(@RHS,[0 tf],Z_o,opts,p);

% unpack Z
th1 = Z(:,1); th2 = Z(:,3); th3 = Z(:,5); th4 = Z(:,7);
th1dot = Z(:,2); th2dot = Z(:,4); th3dot = Z(:,6); th4dot = Z(:,8);

% convert output thetas to cartesian coordinates
% joint positions
x_J1 = p.l1*sin(Z(:,1)); y_J1 = -p.l1*cos(Z(:,1));
x_J2 = x_J1+p.l2*sin(Z(:,3)); y_J2 = y_J1-p.l2*cos(Z(:,3));
x_J3 = x_J2+p.l3*sin(Z(:,5)); y_J3 = y_J2-p.l3*cos(Z(:,5));
x_J4 = x_J3+p.l4*sin(Z(:,7)); y_J4 = y_J3-p.l4*cos(Z(:,7));

% muscle attachment points
x_a11 = p.a11*sin(Z(:,1)); y_a11 = -p.a11*cos(Z(:,1));
x_a12 = x_J1+p.a12*sin(Z(:,3)); y_a12 = y_J1-p.a12*cos(Z(:,3));
x_a22 = x_J1+p.a22*sin(Z(:,3)); y_a22 = y_J1-p.a22*cos(Z(:,3));
x_a23 = x_J2+p.a23*sin(Z(:,5)); y_a23 = y_J2-p.a23*cos(Z(:,5));
x_a33 = x_J2+p.a33*sin(Z(:,5)); y_a33 = y_J2-p.a33*cos(Z(:,5));
x_a34 = x_J3+p.a34*sin(Z(:,7)); y_a34 = y_J3-p.a34*cos(Z(:,7));

% COM of each link, system
x_G1 = p.d1*sin(Z(:,1)); y_G1 = -p.d1*cos(Z(:,1));
x_G2 = x_J1+p.d2*sin(Z(:,3)); y_G2 = y_J1-p.d2*cos(Z(:,3));
x_G3 = x_J2+p.d3*sin(Z(:,5)); y_G3 = y_J2-p.d3*cos(Z(:,5));
x_G4 = x_J3+p.d4*sin(Z(:,7)); y_G4 = y_J3-p.d4*cos(Z(:,7));
x_GM = x_J3+p.l4*sin(Z(:,7)); y_GM = y_J3-p.l4*cos(Z(:,7));

x_COM = (p.m1*x_G1+p.m2*x_G2+p.m3*x_G3+p.m4*x_G4+p.M*x_GM)/mTot;
y_COM = (p.m1*y_G1+p.m2*y_G2+p.m3*y_G3+p.m4*y_G4+p.M*y_GM)/mTot;

%% animate
Time = 0; % Start time of animation
TF = 0.75; % Animation speed, Time viewing = (real time)*TF

% initialize pendulums and axes
pend1 = plot(0,0,'Color','k','LineWidth',2);
pend2 = plot(0,0,'Color','k','LineWidth',2);
pend3 = plot(0,0,'Color','k','LineWidth',2);
pend4 = plot(0,0,'Color','k','LineWidth',2);
mass = rectangle('Position',zeros(1,4));
muscle1 = plot(0,0,'Color','r','LineWidth',3);
muscle2 = plot(0,0,'Color','r','LineWidth',3);
muscle3 = plot(0,0,'Color','r','LineWidth',3);
COM = plot(0,0,'*');
axis([-4.2 4.2 -7.6 0.4]);

tic
% begin animation
while Time < tf/TF
    % interpolate unknown values
    J1 = interp1(t,[x_J1,y_J1],Time*TF);   
    J2 = interp1(t,[x_J2,y_J2],Time*TF);
    J3 = interp1(t,[x_J3,y_J3],Time*TF);
    J4 = interp1(t,[x_J4,y_J4],Time*TF);
    A11 = interp1(t,[x_a11,y_a11],Time*TF);
    A12 = interp1(t,[x_a12,y_a12],Time*TF);
    A22 = interp1(t,[x_a22,y_a22],Time*TF);
    A23 = interp1(t,[x_a23,y_a23],Time*TF);
    A33 = interp1(t,[x_a33,y_a33],Time*TF);
    A34 = interp1(t,[x_a34,y_a34],Time*TF);
    COM_pos = interp1(t,[x_COM,y_COM],Time*TF);

    % define bone geometry
    pend1.XData = [0,J1(1)];
    pend1.YData = [0,J1(2)];
    
    pend2.XData = [J1(1),J2(1)];
    pend2.YData = [J1(2),J2(2)];
    
    pend3.XData = [J2(1),J3(1)];
    pend3.YData = [J2(2),J3(2)];
    
    pend4.XData = [J3(1),J4(1)];
    pend4.YData = [J3(2),J4(2)];
    
    % define body geometry
    mass.Position = [J4(1)-1 J4(2)-1 2 1];
    
    % define muscle geometry
    muscle1.XData = [A11(1),A12(1)];
    muscle1.YData = [A11(2),A12(2)];
    
    muscle2.XData = [A22(1),A23(1)];
    muscle2.YData = [A22(2),A23(2)];
    
    muscle3.XData = [A33(1),A34(1)];
    muscle3.YData = [A33(2),A34(2)];
    
    COM.XData = [COM_pos(1),COM_pos(1)];
    COM.YData = [COM_pos(2),COM_pos(2)];
    
    legend(COM,'Center of mass of system')
    
    % draw
    drawnow;
    Time = toc;
end
%% Plot positions
figure
hold on

plot(x_J1,y_J1)
plot(x_J2,y_J2)
plot(x_J3,y_J3)
plot(x_J4,y_J4)
title('X and Y Position of Each Joint')
xlabel('X Position')
ylabel('Y Position')

figure;
hold on
plot(x_G1,y_G1)
plot(x_G2,y_G2)
plot(x_G3,y_G3)
plot(x_G4,y_G4)
plot(x_GM,y_GM)
plot(x_COM,y_COM)
title('COM position in space')
xlabel('X Position')
ylabel('Y Position')
legend('COM bone 1','COM bone 2','COM bone 3','COM bone 4','COM Body','Total COM')

figure;
hold on
plot(t,sqrt(x_COM.^2+y_COM.^2))
title('X and Y Position of COM vs. Time')
xlabel('Time')
fprintf('Total Movement of Center of Mass = %0.9f m\n',sqrt((x_COM(end)-x_COM(1))^2+(y_COM(end)-y_COM(1))^2))
%% Energy check

% translational kinetic energy of 1st link
trans_vel_x1 = p.d1*cos(th1).*th1dot;
trans_vel_y1 = -p.d1*sin(th1).*th1dot;

trans_vel_squared1 = trans_vel_x1.^2+trans_vel_y1.^2;
trans_kin1 = 0.5*(p.m1)*trans_vel_squared1;

% rotational kinetic energy of first link
rot_kin1 = 0.5*p.I1*th1dot.^2;

% potential energy of first link
potential_p1 = p.m1*p.g*(p.d1-p.d1*cos(th1));

% total energy of first link
total1 = potential_p1+rot_kin1+trans_kin1;

% translational kinetic energy of second link
trans_vel_x2 = p.l1*cos(th1).*th1dot+p.d2*cos(th2).*th2dot;
trans_vel_y2 = -p.l1*sin(th1).*th1dot-p.d2*sin(th2).*th2dot;

trans_vel_squared2 = trans_vel_x2.^2+trans_vel_y2.^2;
trans_kin2 = 0.5*(p.m2)*trans_vel_squared2;

% rotational kinetic energy of second link
rot_kin2 = 0.5*p.I2*th2dot.^2;

% potential energy of second link
potential_p2 = p.m2*p.g*(p.l1+p.d2-(p.l1*cos(th1)+p.d2*cos(th2)));

% total energy of second link
total2 = trans_kin2+rot_kin2+potential_p2;

% translational kinetic energy of third link
trans_vel_x3 = p.l1*cos(th1).*th1dot+p.l2*cos(th2).*th2dot+p.d3*cos(th3).*th3dot;
trans_vel_y3 = -p.l1*sin(th1).*th1dot-p.l2*sin(th2).*th2dot-p.d3*sin(th3).*th3dot;

trans_vel_squared3 = trans_vel_x3.^2+trans_vel_y3.^2;
trans_kin3 = 0.5*(p.m3)*trans_vel_squared3;

% rotational kinetic energy of third link
rot_kin3 = 0.5*p.I3*th3dot.^2;

% potential energy of third link
potential_p3 = p.m3*p.g*(p.l1+p.l2+p.d3-(p.l1*cos(th1)+p.l2*cos(th2)+p.d3*cos(th3)));

% total energy of third link
total3 = trans_kin3+rot_kin3+potential_p3;

% translational kinetic energy of fourth link
trans_vel_x4 = p.l1*cos(th1).*th1dot+p.l2*cos(th2).*th2dot+p.l3*cos(th3).*th3dot+p.d4*cos(th4).*th4dot;
trans_vel_y4 = -p.l1*sin(th1).*th1dot-p.l2*sin(th2).*th2dot-p.l3*sin(th3).*th3dot-p.d4*sin(th4).*th4dot;

trans_vel_squared4 = trans_vel_x4.^2+trans_vel_y4.^2;
trans_kin4 = 0.5*(p.m4)*trans_vel_squared4;

% rotational kinetic energy of fourth link
rot_kin4 = 0.5*p.I4*th4dot.^2;

% potential energy of fourth link
potential_p4 = p.m4*p.g*(p.l1+p.l2+p.l3+p.d4-(p.l1*cos(th1)+p.l2*cos(th2)+p.l3*cos(th3)+p.d4*cos(th4)));

% total energy of fourth link
total4 = trans_kin4+rot_kin4+potential_p4;
kinetic_4 = trans_kin4+rot_kin4;

% translational kinetic energy of M
trans_vel_xM = p.l1*cos(th1).*th1dot+p.l2*cos(th2).*th2dot+p.l3*cos(th3).*th3dot+p.l4*cos(th4).*th4dot;
trans_vel_yM = -p.l1*sin(th1).*th1dot-p.l2*sin(th2).*th2dot-p.l3*sin(th3).*th3dot-p.l4*sin(th4).*th4dot;

trans_vel_squaredM = trans_vel_xM.^2+trans_vel_yM.^2;
trans_kinM = 0.5*(p.M)*trans_vel_squaredM;

% potential energy of M link
potential_pM = p.M*p.g*(p.l1+p.l2+p.l3+p.l4-(p.l1*cos(th1)+p.l2*cos(th2)+p.l3*cos(th3)+p.l4*cos(th4)));

% total energy of M link
totalM = trans_kinM+potential_pM;

% total system energy
TOTAL = total1+total2+total3+total4+totalM;

% plot
figure
hold on
plot(t,total1)
plot(t,total2)
plot(t,total3)
plot(t,total4)
plot(t,totalM)
plot(t,TOTAL)
title('Energy vs. Time')
xlabel('Time')
ylabel('Energy')
legend('Link 1','Link 2','Link 3','Link 4','Mass M','Total System Energy')
end
function [Zdot] = RHS(t,Z,p)

% unpack Z vector
th1 = Z(1); th1dot = Z(2); th2 = Z(3); th2dot = Z(4); th3 = Z(5); th3dot = Z(6);
th4 = Z(7); th4dot = Z(8);

% Send to MATLAB built acceleration solver (min for minimal coordinates,
% DAE for DAE solution, Lag for lagrange)
[A,B] = IA_accel(p.I1,p.I2,p.I3,p.I4,p.M,p.M1,p.M2,p.M3,...
    p.a11,p.a12,p.a22,p.a23,p.a33,p.a34,p.d1,p.d2,p.d3,p.d4,p.g,...
    p.l1,p.l2,p.l3,p.l4,p.m1,p.m2,p.m3,p.m4,...
    th1,th2,th3,th4,th1dot,th2dot,th3dot,th4dot);

x = A\B;

% Pack Zdot vector
Zdot = [th1dot x(11) th2dot x(12) th3dot x(13) th4dot x(14)]';
end