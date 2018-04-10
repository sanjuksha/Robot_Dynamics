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

for i=1:6
    % For symbolic representations from T0 to T6
    M = dhparam2matrix(d(i),t(i), a(i), b(i))
    H = H * M;
end
