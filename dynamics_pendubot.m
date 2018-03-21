%Created on: March, 2018
%Author(s): Pradheep Krishna Muthukrishnan Padmanabhan
%Copyright (c) [2018]
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.

% This file contains the dynamic model of the pendubot. So based on this
% file a swing up control using feedback linearization will be developed
%This otherwise also called as computed torque technique.

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
%j =jacobian([p_q_1,p_y1,p_q_2,p_y2],[q_1,q_2]);
%v = j*q;
v1x = -l_1*q_1d*sin(q_1)
v1y = l_1*q_1d*cos(q_1)
v2x = v1x-l_2*sin(q_1+q_2)*(q_1d+q_2d);
v2y = v1y+l_2*cos(q_1+q_2)*(q_1d+q_2d);
%P and V are the velocities that have been obtained

%Step 3, Now calculating the Kinetic energy and the potential energy based
%on the lagrangian mechanics
KE_1 = 0.5*m1*(v1x+v1y)^2;
KE_2 = 0.5*m2*(v2x+v2y)^2;

KE = KE_1 + KE_2;
PE = m1*g*p_y1 +m2*g*p_y2;
KE = simplify(KE);
PE = simplify(PE);
L = KE-PE; 
matlabFunction(L, 'file', 'Lagrange');

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

%In the below step, We establish a linear relationship between output and
%U.
syms kd kp q_1ddr q_1dr q_1r %A PD control is being introduced in order for the tracking control 
%which is based on the outer loop control, See schematic for better
%understanding

%In order to track, we setup a reference trajectory! 
v_1 = q_1ddr+kd*(q_1dr-q_1d) + kp*(q_1r - q_1);
tau_1 = simplify( ddtpKEpq_1d - pKEpq_1 + pPEpq_1-p_q_1)==0;
tau_2 = simplify( ddtpKEpq_2d - pKEpq_2 + pPEpq_2 - p_q_2)==0;

%For the state space equation
q_1dds = isolate(tau_1,q_1dd);
q_2dds = isolate(tau_1,q_2dd);
q_1dds01 = isolate(tau_2,q_1dd);
q_2dds01 = isolate(tau_2,q_2dd);
%Hard coding the expression for state space transformation
q_2ddss1 = (l_1^2*m1*q_1dd*sin(2*q_1) - l_1^2*m2*q_1dd - l_2^2*m2*q_1dd - g*l_2*m2*cos(q_1 + q_2) - l_1^2*m1*q_1dd + l_1^2*m2*q_1dd*sin(2*q_1) - g*l_1*m1*cos(q_1) - g*l_1*m2*cos(q_1) + l_2^2*m2*q_1dd*sin(2*q_1 + 2*q_2) + l_1^2*m1*q_1d^2*cos(2*q_1) + l_1^2*m2*q_1d^2*cos(2*q_1) + l_2^2*m2*q_1d^2*cos(2*q_1 + 2*q_2) + l_2^2*m2*q_2d^2*cos(2*q_1 + 2*q_2) + 2*l_1*l_2*m2*q_1dd*sin(2*q_1 + q_2) + l_1*l_2*m2*q_2d^2*sin(q_2) + 2*l_1*l_2*m2*q_1d^2*cos(2*q_1 + q_2) + l_1*l_2*m2*q_2d^2*cos(2*q_1 + q_2) + 2*l_2^2*m2*q_1d*q_2d*cos(2*q_1 + 2*q_2) - 2*l_1*l_2*m2*q_1dd*cos(q_2) + 2*l_1*l_2*m2*q_1d*q_2d*sin(q_2) + 2*l_1*l_2*m2*q_1d*q_2d*cos(2*q_1 + q_2))/(l_2^2*m2 - l_2^2*m2*sin(2*q_1 + 2*q_2) + l_1*l_2*m2*cos(q_2) - l_1*l_2*m2*sin(2*q_1 + q_2));
q_1ddss1 = (l_2^2*m2*q_2dd*sin(2*q_1 + 2*q_2) - g*l_2*m2*cos(q_1 + q_2) - g*l_1*m1*cos(q_1) - g*l_1*m2*cos(q_1) - l_2^2*m2*q_2dd + l_1^2*m1*q_1d^2*cos(2*q_1) + l_1^2*m2*q_1d^2*cos(2*q_1) + l_2^2*m2*q_1d^2*cos(2*q_1 + 2*q_2) + l_2^2*m2*q_2d^2*cos(2*q_1 + 2*q_2) + l_1*l_2*m2*q_2dd*sin(2*q_1 + q_2) + l_1*l_2*m2*q_2d^2*sin(q_2) + 2*l_1*l_2*m2*q_1d^2*cos(2*q_1 + q_2) + l_1*l_2*m2*q_2d^2*cos(2*q_1 + q_2) + 2*l_2^2*m2*q_1d*q_2d*cos(2*q_1 + 2*q_2) - l_1*l_2*m2*q_2dd*cos(q_2) + 2*l_1*l_2*m2*q_1d*q_2d*sin(q_2) + 2*l_1*l_2*m2*q_1d*q_2d*cos(2*q_1 + q_2))/(l_1^2*m1 + l_1^2*m2 + l_2^2*m2 - l_2^2*m2*sin(2*q_1 + 2*q_2) - l_1^2*m1*sin(2*q_1) - l_1^2*m2*sin(2*q_1) + 2*l_1*l_2*m2*cos(q_2) - 2*l_1*l_2*m2*sin(2*q_1 + q_2));
%Tau_2 should be zero 
q_1ddss21 = (l_2*cos(q_1 + q_2) + l_1*cos(q_1) - l_2^2*m2*q_2dd - g*l_2*m2*cos(q_1 + q_2) + l_2^2*m2*q_2dd*sin(2*q_1 + 2*q_2) + l_2^2*m2*q_1d^2*cos(2*q_1 + 2*q_2) + l_2^2*m2*q_2d^2*cos(2*q_1 + 2*q_2) - l_1*l_2*m2*q_1d^2*sin(q_2) + l_1*l_2*m2*q_1d^2*cos(2*q_1 + q_2) + 2*l_2^2*m2*q_1d*q_2d*cos(2*q_1 + 2*q_2))/(l_2^2*m2 - l_2^2*m2*sin(2*q_1 + 2*q_2) + l_1*l_2*m2*cos(q_2) - l_1*l_2*m2*sin(2*q_1 + q_2));
q_2ddss21 = (l_2*cos(q_1 + q_2) + l_1*cos(q_1) - l_2^2*m2*q_1dd - g*l_2*m2*cos(q_1 + q_2) + l_2^2*m2*q_1dd*sin(2*q_1 + 2*q_2) + l_2^2*m2*q_1d^2*cos(2*q_1 + 2*q_2) + l_2^2*m2*q_2d^2*cos(2*q_1 + 2*q_2) + l_1*l_2*m2*q_1dd*sin(2*q_1 + q_2) - l_1*l_2*m2*q_1d^2*sin(q_2) + l_1*l_2*m2*q_1d^2*cos(2*q_1 + q_2) + 2*l_2^2*m2*q_1d*q_2d*cos(2*q_1 + 2*q_2) - l_1*l_2*m2*q_1dd*cos(q_2))/(l_2^2*m2 - l_2^2*m2*sin(2*q_1 + 2*q_2));
q_1ddss2 = subs(q_1ddss21,tau_2,0);
q_2ddss2 = subs(q_2ddss21,tau_2,0);

%Continuing the model
eqq_2 = simplify( ddtpKEpq_2d - pKEpq_2 + pPEpq_2)==0; 
qSol = isolate(eqq_2,q_2dd);
q_2dd1 = subs(eqq_2, lhs(qSol), rhs(qSol));
q_2dd2 = -l_2*m2*(l_2*q_1d^2*cos(2*q_1 + 2*q_2) - g*cos(q_1 + q_2) - l_2*q_1dd + l_2*q_2d^2*cos(2*q_1 + 2*q_2) - l_1*q_1dd*cos(q_2) + l_1*q_1dd*sin(2*q_1 + q_2) - l_1*q_1d^2*sin(q_2) - (l_2*(l_2*q_1d^2*cos(2*q_1 + 2*q_2) - g*cos(q_1 + q_2) - l_2*q_1dd + l_2*q_2d^2*cos(2*q_1 + 2*q_2) - l_1*q_1dd*cos(q_2) + l_1*q_1dd*sin(2*q_1 + q_2) - l_1*q_1d^2*sin(q_2) + l_1*q_1d^2*cos(2*q_1 + q_2) + l_2*q_1dd*sin(2*q_1 + 2*q_2) + 2*l_2*q_1d*q_2d*cos(2*q_1 + 2*q_2)))/(l_2 - l_2*sin(2*q_1 + 2*q_2)) + l_1*q_1d^2*cos(2*q_1 + q_2) + l_2*q_1dd*sin(2*q_1 + 2*q_2) + (l_2*sin(2*q_1 + 2*q_2)*(l_2*q_1d^2*cos(2*q_1 + 2*q_2) - g*cos(q_1 + q_2) - l_2*q_1dd + l_2*q_2d^2*cos(2*q_1 + 2*q_2) - l_1*q_1dd*cos(q_2) + l_1*q_1dd*sin(2*q_1 + q_2) - l_1*q_1d^2*sin(q_2) + l_1*q_1d^2*cos(2*q_1 + q_2) + l_2*q_1dd*sin(2*q_1 + 2*q_2) + 2*l_2*q_1d*q_2d*cos(2*q_1 + 2*q_2)))/(l_2 - l_2*sin(2*q_1 + 2*q_2)) + 2*l_2*q_1d*q_2d*cos(2*q_1 + 2*q_2));
tau_1 = subs(tau_1,q_2dd,q_2dd2);

tau_1 = subs(tau_1,q_1dd,v_1)
matlabFunction(tau_1, 'file', 'ControlTorque1');
%Internal dynamics of the system and their stability

%tau_2 should be zero
