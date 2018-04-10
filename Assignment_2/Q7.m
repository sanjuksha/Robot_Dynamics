clc; 
% Input D-H parameters
% alpha or link twist = b
b =  [-pi/2 0 -pi/2 pi/2 -pi/2 0];

% Link length = a
a =[0 270 70 0 0 0];
% Input Joint angles theta as t
t =  [pi/6 pi/3-pi/2 pi/6 -pi/5 pi/12 pi/5+pi];

% Link offset = d
d = [290 0 0 302 0 72];
% Joint angle = theta = t(i)

T01 = dhparam2matrix(d(1),t(1), a(1),b(1));

T12 = dhparam2matrix(d(2),t(2),a(2),b(2));

T23 = dhparam2matrix(d(3),t(3),a(3),b(3));

T34 = dhparam2matrix(d(4),t(4),a(4),b(4));

T45 = dhparam2matrix(d(5),t(5),a(5),b(5));

T56 = dhparam2matrix(d(6),t(6),a(6),b(6));

T06 = T01*T12*T23*T34*T45*T56

T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;

x = [T01(1,4) T02(1,4) T03(1,4) T04(1,4) T05(1,4) T06(1,4)];
y = [T01(2,4) t02(2,4) t03(2,4) t04(2,4) T05(2,4) T06(2,4)];
z = [T01(3,4) T02(3,4) T03(3,4) T04(3,4) T05(3,4) T06(3,4)]; 

px = T06*[70;0;0;1];
py = T06*[0;70;0;1];
pz = T06*[0;0;70;1];

px1 = [T06(1,4),px(1,1)] ;
py1 = [T06(2,4),px(2,1)]; 
pz1 = [T06(3,4),px(3,1)] ;

px2 = [T06(1,4),py(1,1)] ;
py2 = [T06(2,4),py(2,1)] ;
pz2 = [T06(3,4),py(3,1)] ;

px3 = [T06(1,4),pz(1,1)] ;
py3 = [T06(2,4),pz(2,1)] ;
pz3 = [T06(3,4),pz(3,1)] ;


plot3(x,y,z,'-o');
hold on 
%plot3(ptx1,pty1,ptz1,'r','LineWidth',3)
%hold on 
%plot3(ptx2,pty2,ptz2,'g','LineWidth',3)
%hold on 
%plot3(ptx3,pty3,ptz3,'b','LineWidth',3)
title("Forward Kinematics")
xlabel("x(mm)")
ylabel("y(mm)")
zlabel("z(mm)")
plot3(x,y,z, 'y','Linewidth',3);



