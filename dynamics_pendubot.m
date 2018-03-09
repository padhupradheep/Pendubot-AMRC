% This file contains the dynamic model of the pendubot. So based on this
% file a swing up control using feedback linearization will be developed

% Step 1. Symbolic representation of position and velocity for lagrangian mechanics

clear;close all;clc;
syms m1 m2 l_1 l_2 q_1 q_2 q_1d q_2d; %Pendubot's mass of the links, length of the links and the joint angles
syms tau_1 tau_2; %Torque
syms g; 
syms q_1d q_1dd q_2d q_2dd;
syms v_1;

% Step 2. Now computing the position and velocities of the masses in the
% system:
y = [q_1;q_2];
p_q_1 = l_1*cos(q_1);  
p_y1 = l_1*sin(q_1);

p_q_2 = p_q_1+l_2*cos(q_1+q_2);
p_y2 = p_y1+l_2*sin(q_1+q_2);

p = [p_q_1 p_y1;
    p_q_2 p_y2];
q = [q_1d;q_2d];
j =jacobian([p_q_1,p_y1,p_q_2,p_y2],[q_1,q_2]);
v = j*q;


%P and V are the velocities that have been obtained

%Step 3, Now calculating the Kinetic energy and the potential energy based
%on the lagrangian mechanics
KE_1 = 0.5*m1*(v(1,1)+v(2,1))^2;
KE_2 = 0.5*m2*(v(3,1)+v(4,1))^2;

KE = KE_1 + KE_2;
PE = m1*g*p_y1 +m2*g*p_y2;
KE = simplify(KE);
PE = simplify(PE);
L = KE-PE; 
matlabFunction(L, 'file', 'Lagrange');
%Deriving the equaations of motion from Lagrangian: 

pKEpq_1d = diff(KE,q_1d);
ddtpKEpq_1d = diff(pKEpq_1d,q_1)*q_1d+ ...
             diff(pKEpq_1d,q_1d)*q_1dd+ ...
             diff(pKEpq_1d,q_2)*q_2d + ...
             diff(pKEpq_1d,q_2d)*q_2dd;
pKEpq_1 = diff(KE,q_1);
pPEpq_1 = diff(PE,q_1);

%Using chain rule in order to derive the lagrange mechanics
pKEpq_2d = diff(KE,q_2d);
ddtpKEpq_2d = diff(pKEpq_2d,q_1)*q_1d+ ...
             diff(pKEpq_2d,q_1d)*q_1dd+ ...
             diff(pKEpq_2d,q_2)*q_2d + ...
             diff(pKEpq_2d,q_2d)*q_2dd;
pKEpq_2 = diff(KE,q_2);
pPEpq_2 = diff(PE,q_2);   

syms kd kp 
v_1 = kd*-q_1d + kp*(q_1d - q_1);
tau_1 = simplify( ddtpKEpq_1d - pKEpq_1 + pPEpq_1);
tau_1 = subs(tau_1,q_1dd,v_1)
matlabFunction(tau_1, 'file', 'ControlTorque1');
eqq_2 = simplify( ddtpKEpq_2d - pKEpq_2 + pPEpq_2) == 0 
%tau_2 should be zero

q_2dd = isolate(eqq_2,q_2dd)
q_2dd
