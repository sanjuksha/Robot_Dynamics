
%% Question 1- Forward Kinematics to find the tip position
syms q1 q2 q3 L1 L2 L3 mL g m1 m2 m3 h1 h2 h3 hL q_d1 q_d2 q_d3
clc;

% Rotation
R0 = rotate(q1);
R1 = rotate(q2);
R2 = rotate(q3);

% Translation

T01 = translate(L1);
T12 = translate(L2);
T23 = translate(L3);

% Final Transformation Matrix
R01 =R0*T01;
R02 = R01*R1*T12;
R03 = R01*R02*R2*T23;

T03 = simplify(R0*T01*R1*T12*R2*T23)

% End Effector Position
EEP = T03(1:3,4)
EEO = EulerAngle(T03)
%% Question 2 - Velocity Kinematics
J = [diff(T03(1,4),q1) diff(T03(1,4),q2) diff(T03(1,4),q3); 
    diff(T03(2,4),q1) diff(T03(2,4),q2) diff(T03(2,4),q3);
    diff(T03(3,4),q1) diff(T03(3,4),q2) diff(T03(3,4),q3);
    R01(1,1) R02(1,1) R03(1,1) ;
    R01(2,1) R02(2,1) R03(2,1);
    R01(3,1) R02(3,1) R03(3,1)];
%% Question 3a
% Tip Position
EEP_val = double(subs(EEP,[q1,q2,q3,L1,L2,L3],[deg2rad(30),deg2rad(20),deg2rad(-25),0.6,0.4,0.15]))
disp('Units: m ')
% Tip Orientation
EEO_val = double(subs(EEO,[q1,q2,q3,L1,L2,L3],[deg2rad(30),deg2rad(20),deg2rad(-25),0.6,0.4,0.15]))
disp('Units : deg')

%% Question 3b
Jacobian = double(subs(J,[q1,q2,q3,L1,L2,L3],[deg2rad(30),deg2rad(20),deg2rad(-25),0.6,0.4,0.15]))

%% Question 4
J_T = transpose(J);
F = [0;0;-mL*g;0;0;0];
% Symbolic Torque
t = simplify(J_T*F)
%Numerical Torque
Torque = double(subs(t,[q1,q2,q3,L1,L2,L3,mL,g],[deg2rad(30),deg2rad(20),deg2rad(-25),0.6,0.4,0.15,1.5,9.81]))
%% 
syms  q1(t) q2(t) q3(t) y1(t) y2(t) y3(t) yL(t) z1(t) z2(t) z3(t) zL(t)
q_d1 = diff(q1);
q_d2 = diff(q2);
q_d3 = diff(q3);
q_dot = [q_d1;q_d2;q_d3];

x_dot = J*q_dot;

%% Question 3c

velocity = double(Jacobian * subs(q_dot,[q_d1,q_d2,q_d3],[deg2rad(45),deg2rad(45),deg2rad(45)]));
disp('Units : deg/sec')
%% Question 5
% To find Lagrange's equation 
% Potential Energy of Each Link

y1 = (L1/2)* cos(q1);
z1 = (L1/2)* sin(q1);
P1 = potentialEnergy(m1,z1);

y2 = L1*cos(q1) + (L2/2)*cos(q1+q2);
z2 = L1*sin(q1) + (L2/2)*sin(q1+q2);
P2 = potentialEnergy(m2,z2);

y3 = L1*cos(q1) + L2*cos(q1+q2) + (L3/2)*cos(q1+q2+q3);
z3 = L1*sin(q1) + L2*sin(q1+q2) + (L3/2)*sin(q1+q2+q3);
P3 = potentialEnergy(m3,z3);

yL = L1*cos(q1) + L2*cos(q1+q2) + L3*cos(q1+q2+q3);
zL = L1*sin(q1) + L2*sin(q1+q2) + L3*sin(q1+q2+q3);
PL = potentialEnergy(mL,zL);

P = P1 + P2 + P3 + PL


% Kinetic energy of each link 
m_1_P = [0; y1;z1];
v1 = diff(m_1_P);
K1 = kineticEnergy(m1,v1);

m_2_P = [0; y2;z2];
v2 = diff(m_2_P);
K2 = kineticEnergy(m2,v2);


m_3_P = [0; y3;z3];
v3 = diff(m_3_P);
K3 = kineticEnergy(m3,v3);


m_3_P = [0; y3;z3];
vL = diff(m_3_P);
KL = kineticEnergy(mL,vL);
K = K1 + K2 + K3 + KL

% Lagrangian Equation

syms q1_dot q2_dot q3_dot q1n q2n q3n t q1_ddot q2_ddot q3_ddot
L = K - P;
L_val = simplify(subs(L,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8]))
L_n = simplify(subs(L,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)],[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n]));

% Torque Equations

% FOR TAU 1
term2 = diff(L_n,q1n);
term1 = diff(L_n,q1_dot);
L_r = subs(term1,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm2 = subs(term2,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm1 = diff(L_r);
final1 = fterm1 - fterm2;
tau = double(simplify(subs(final1,[L1,L2,L3,m1,m2,m3,mL,q1(t),q2(t),q3(t),g],[0.6,0.4,0.15,3,2,0.5,2,deg2rad(30),deg2rad(20),deg2rad(-25),9.8])));

% FOR TAU 2
term2s = diff(L_n,q2n);
term1s = diff(L_n,q2_dot);
L_r2 = subs(term1,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm2s = subs(term2s,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm1s = diff(L_r2);
final2 = fterm1s - fterm2s;
tau2 = double(simplify(subs(final2,[L1,L2,L3,m1,m2,m3,mL,q1(t),q2(t),q3(t),g],[0.6,0.4,0.15,3,2,0.5,2,deg2rad(30),deg2rad(20),deg2rad(-25),9.8])));

% FOR TAU 3
term2t = diff(L_n,q3n);
term1t = diff(L_n,q3_dot);
L_r3 = subs(term1t,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm2t = subs(term2t,[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n],[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t)]);
fterm1t = diff(L_r3);
final3 = fterm1t - fterm2t;
tau3 = double(simplify(subs(final3,[L1,L2,L3,m1,m2,m3,mL,q1(t),q2(t),q3(t),g],[0.6,0.4,0.15,3,2,0.5,2,deg2rad(30),deg2rad(20),deg2rad(-25),9.8])));

% Simplifying symbolic TAU
tau_sym1 = subs(final1,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t),diff(q1(t),t,t),diff(q2(t),t,t),diff(q3(t),t,t)],[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n, q1_ddot, q2_ddot,q3_ddot]);
tau_sym2 = subs(final2,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t),diff(q1(t),t,t),diff(q2(t),t,t),diff(q3(t),t,t)],[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n, q1_ddot, q2_ddot,q3_ddot]);
tau_sym3 = subs(final3,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),q1(t),q2(t),q3(t),diff(q1(t),t,t),diff(q2(t),t,t),diff(q3(t),t,t)],[q1_dot, q2_dot, q3_dot, q1n, q2n, q3n, q1_ddot, q2_ddot,q3_ddot]);

Tau_sym_d =[tau_sym1; tau_sym2; tau_sym3];
Tau_sym_d = Tau_sym_d(t);

% Value of TAU according to give configuration
Torque = [tau ; tau2; tau3];

% Gravity Matrix
g1 = simplify(subs(tau_sym1,[q1_ddot,q1_dot,q2_ddot,q2_dot,q3_ddot,q3_dot],[0,0,0,0,0,0]));
g2 = simplify(subs(tau_sym2,[q1_ddot,q1_dot,q2_ddot,q2_dot,q3_ddot,q3_dot],[0,0,0,0,0,0]));
g3 = simplify(subs(tau_sym3,[q1_ddot,q1_dot,q2_ddot,q2_dot,q3_ddot,q3_dot],[0,0,0,0,0,0]));
G = [g1;g2;g3];

% M Matrix

new_tau11 =  simplify(subs(Tau_sym_d(1),q1_ddot,0));
tau_m11 = simplify(Tau_sym_d(1) - new_tau11);
t11 = tau_m11/q1_ddot;

new_tau12 =  simplify(subs(Tau_sym_d(1),q2_ddot,0));
tau_m12 = simplify(Tau_sym_d(1) - new_tau12);
t12 = tau_m12/q2_ddot;

new_tau13 =  simplify(subs(Tau_sym_d(1),q3_ddot,0));
tau_m13 = simplify(Tau_sym_d(1) - new_tau13);
t13 = tau_m13/q3_ddot;

new_tau21 =  simplify(subs(Tau_sym_d(2),q1_ddot,0));
tau_m21 = simplify(Tau_sym_d(2) - new_tau21);
t21 = tau_m21/q1_ddot;

new_tau22 =  simplify(subs(Tau_sym_d(2),q2_ddot,0));
tau_m22 = simplify(Tau_sym_d(2) - new_tau22);
t22 = tau_m22/q2_ddot;

new_tau23 =  simplify(subs(Tau_sym_d(2),q3_ddot,0));
tau_m23 = simplify(Tau_sym_d(2) - new_tau23);
t23 = tau_m23/q3_ddot;

new_tau31 =  simplify(subs(Tau_sym_d(3),q1_ddot,0));
tau_m31 = simplify(Tau_sym_d(3) - new_tau31);
t31 = tau_m31/q1_ddot;

new_tau32 =  simplify(subs(Tau_sym_d(3),q2_ddot,0));
tau_m32 = simplify(Tau_sym_d(3) - new_tau32);
t32 = tau_m32/q2_ddot;

new_tau33 =  simplify(subs(Tau_sym_d(3),q3_ddot,0));
tau_m33 = simplify(Tau_sym_d(3) - new_tau33);
t33 = tau_m33/q3_ddot;

M_matrix =[t11 t12 t13; t21 t22 t23; t31 t32 t33]
M_value = vpa(simplify(subs(M_matrix,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])),2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C MATRIX %%%%%%%%%%%%%%%

q_ddot = [q1_ddot; q2_ddot;q3_ddot];
q_dot = [q1_dot; q2_dot;q3_dot];
% C MATRIX
M_qd = M_matrix*q_ddot;
C_m = simplify((Tau_sym_d - (M_qd + G)),50);
C_m =C_m(t);
C_m11 = C_m(1)/q1_dot;
C_m12 = C_m(1)/q2_dot;
C_m13 = C_m(1)/q3_dot;

C_m21 = C_m(2)/q1_dot;
C_m22 = C_m(2)/q2_dot;
C_m23 = C_m(2)/q3_dot;

C_m31 = C_m(3)/q1_dot;
C_m32 = C_m(3)/q2_dot;
C_m33 = C_m(3)/q3_dot;

C_matrix = [C_m11 C_m12 C_m13; C_m21 C_m22 C_m23; C_m31 C_m32 C_m33];
C_value =simplify(subs(C_matrix,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8]))

%% STANDARD FORM

eqn = vpa(simplify(Tau_sym_d == M_value*q_ddot + C_value*q_dot + G),2)

function [K]= kineticEnergy(m,v)
    K = 0.5 * m * transpose(v)*v;
end

function [P]= potentialEnergy(m,h)
    global g
    P = m*g*h;
end

function [R] =rotate(th)
    R =[1 0 0  0; 0 cos(th) -sin(th) 0; 0 sin(th) cos(th) 0; 0 0 0 1 ];
end

function [T]= translate(L)
    T = [1 0 0 0; 0 1 0 L; 0 0 1 0; 0 0 0 1];
end
function rot_angle = EulerAngle(Tfinal)
% Convert rotation matrix to Euler angles
singularity_check = sqrt(Tfinal(1,1)^2 + Tfinal(2,1)^2);
if singularity_check < 10^-6
    PHI = atan2(Tfinal(2,3), Tfinal(2,2)); % corresponds to RX in Robot Studio
    THETA = atan2(-Tfinal(3,1), sqrt(Tfinal(1,1)^2+Tfinal(2,1)^2)); % corresponds to RY in Robot Studio
    PSI = 0; % corresponds to RZ in Robot Studio
else
    THETA = atan2(-Tfinal(3,1),sqrt(Tfinal(1,1)^2+Tfinal(2,1)^2)); % corresponds to RY in Robot Studio
    PHI = atan2(Tfinal(3,2),Tfinal(3,3)); % corresponds to RX in Robot Studio
    PSI = atan2(Tfinal(2,1),Tfinal(1,1)); % corresponds to RZ in Robot Studio
end

rot_angle = [PHI,THETA,PSI]*180/pi;
end