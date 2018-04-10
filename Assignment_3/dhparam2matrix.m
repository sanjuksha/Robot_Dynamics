function M = dhparam2matrix(l,t, a, alpha)
old=digits(3);
 M = [cos(t) (-sin(t).*cos(alpha)) (sin(t).*sin(alpha)) (a.*cos(t));
     sin(t) cos(t).*cos(alpha) (-cos(t).*sin(alpha)) a.*sin(t);
     0 (sin(alpha)) cos(alpha) l; 0 0 0 1];
end
% Matlab function that takes set of D-H parameters as nx1 vector and
% returns 4x4 transformation matrix
