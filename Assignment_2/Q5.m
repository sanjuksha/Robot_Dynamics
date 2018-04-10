clc;
syms q1 q2 q3 q4 q5 q6
% q1, q2,q3,q4,q5,q6 = 0
% Input D-H parameters
% alpha or link twist = b
b = [-pi/2 0 -pi/2 pi/2 -pi/2 0];

% Link length = a
a = [0 270 70 0 0 0];
% Input Joint angles theta as t
t =[0 (-pi/2) 0 0 0 pi];
% Link offset = d
d = [290 0 0 302 0 72];
H = eye(4);
% Origin = o 
o = [0;0;0;1]


%From question 4
for i=1:6
    % For symbolic representations from T0 to T6 
    % plugged in [0,0,0,0,0,0]
    M = dhparam2matrix(d(i),t(i), a(i), b(i))
    H = H * M;
end
% Rotation matrix 
Rotation = H(1:3,1:3)
% x,y,z co-ordinates
EE_Position = H*o

% the x ,y,z co-ordinates matches the one from Robostudio
