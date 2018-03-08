% This file contains the dynamic model of the pendubot. So based on this
% file a swing up control using feedback linearization will be developed

% Step 1. Symbolic representation of position and velocity for lagrangian mechanics

clear;close all;clc;
syms m1 m2 l_1 l_2 q_1 q_2 q_1d q_2d; %Pendubot's mass of the links, length of the links and the joint angles
syms tau_1 tau_2; %Torque
syms g; 

% Step 2. Now computing the position and velocities of the masses in the
% system:

p_x1 = l_1*cos(q_1);  
p_y1 = l_1*sin(q_1);

p_x2 = p_x1+l_2*cos(q_1+q_2);
p_y2 = p_y1+l_2*sin(q_1+q_2);

p = [p_x1 p_y1;
    p_x2 p_y2];
q = [q_1d;q_2d];
j =jacobian([p_x1,p_y1,p_x2,p_y2],[q_1,q_2]);
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
L = KE-PE

