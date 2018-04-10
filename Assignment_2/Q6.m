clc; 
% Input D-H parameters
% alpha or link twist = b
b = [-pi/2 0 -pi/2 pi/2 -pi/2 0];

% Link length = a
a = [0 270 70 0 0 0];
% Input Joint angles theta as t
t = [-pi/4 (pi/6-pi/2) -pi/6 -pi/6 -pi/4 (pi+pi)];

% Link offset = d
d = [290 0 0 302 0 72];
% Joint angle = theta = t(i)
% H = identity matrix to start multiplication
H = eye(4);
o = [0;0;0;1];


for i=1:6
    M = dhparam2matrix(d(i),t(i), a(i), b(i));
    H = H * M;
end

display(H)
EE_Position = H*o

% Cartesian Position 
c_pos = H(1:3,4)
% Approach Vector
AV = H(1:3,3)

