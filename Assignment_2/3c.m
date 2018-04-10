syms q1 q2 q3 q4 q5 q6
clc; 
% Input D-H parameters
% alpha or link twist = b
b = [-pi/2 0 -pi/2 pi/2 -pi/2 0];

% Link length = a
a = [0 270 70 0 0 0];
% Input Joint angles theta as t
t =[q1 (q2-pi/2) q3 q4 q5 (q6+pi)];
% Link offset = d
d = [290 0 0 302 0 72];

T01 = dhparam2matrix(t(1),d(1),a(1),b(1));

T12 = dhparam2matrix(t(2),d(2),a(2),b(2));

T23 = dhparam2matrix((3),d(3),a(3),b(3));

T34 = dhparam2matrix(t(4),d(4),a(4),b(4));

T45 = dhparam2matrix(t(5),d(5),a(5),b(5));

T56 = dhparam2matrix(t(6),d(6),a(6),b(6));


% Intermediate Transformations
T02 = T01*T12
T03 = T02*T23
T04 = T03*T34
T05 = T04*T45
T06 = T05*T56