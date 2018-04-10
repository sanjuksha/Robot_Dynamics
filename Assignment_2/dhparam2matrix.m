function M = dhparam2matrix(d,t, a, b)
 M = [cos(t) (-sin(t).*cos(b)) (sin(t).*sin(b)) (a.*cos(t));
     sin(t) cos(t).*cos(b) (-cos(t).*sin(b)) a.*sin(t);
     0 (sin(b)) cos(b) d; 0 0 0 1];
end
% Matlab function that takes set of D-H parameters as nx1 vector and
% returns 4x4 transformation matrix