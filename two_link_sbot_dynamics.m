syms q1 q2 q3 dq1 dq2 dq3 mt ml mw Iw l d lt r g real
syms hd xd q1d q2d q3d kpx kph kdx kdh kp1 kp2 kp3 kd1 kd2 kd3 psid real

% relative coords
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
s = [q;dq];

% absolute transformation
T = [1 0 0;
     1 1 0;
     1 1 1];

th = T*q - [0; 0; pi];
dth = T*dq;

a = th(2)-pi; % angle from wheel vertical to leg

% world positions
p1 = [r*th(3); r];                  % wheel center pos
p2 = p1 + [l*sin(a); l*cos(a)];     % leg-torso joint pos
p3 = p2 + [lt*sin(q1); lt*cos(q1)]; % torso pos
p4 = p1 + [d*sin(a); d*cos(a)];     % leg mass pos

pe = p1 + r*[sin(th(3)); cos(th(3))]; % encoder pos

% world velocities
dp1 = [r*dth(3); 0];                          % wheel center
dp2 = dp1 + dth(2)*[l*cos(a); -l*sin(a)];     % leg-torso joint
dp3 = dp2 + dq1*[lt*cos(q1); -lt*sin(q1)];    % torso
dp4 = dp1 + dth(2)*[d*cos(a); -d*sin(a)];     % leg mass

disp('Kinematics Complete.');

% Lagrangian Dynamics
KE = .5*dp1'*mw*dp1 + .5*Iw*dth(3)^2 + .5*dp4'*ml*dp4 + .5*dp3'*mt*dp3; 
PE = mw*g*r + ml*g*p4(2) + mt*g*p3(2);
qAct = [q2;q3]; % actuated angles

[D, C, G, B] = LagrangianDynamics(KE, PE, q, dq, qAct);

% Inverse dynamics
f_s = simplify([dq; D\(-C*dq-G)]);
g_s = simplify([zeros(size(q,1),size(qAct,1)); D\B]);

disp('Dynamics Complete.');

%% Output Linearization

pCom = (ml*p4+mt*p3)/(ml+mt);

psi = atan(pCom(2)/pCom(1));

y = [psi-psid]; % control output

Lfy = simplify(jacobian(y,s)*f_s);  % = dydt
Lgy = simplify(jacobian(y,s)*g_s);  % should be 0
LgLfy = simplify(jacobian(Lfy,s)*g_s);
Lf2y = simplify(jacobian(Lfy,s)*f_s);

disp('Lie Derivatives Complete.');

Kp = diag([kp1 kp2 kp3]); Kd = diag([kd1 kd2 kd3]); % PD controller gains

u = simplify(LgLfy\(-Lf2y -Kp*y -Kd*Lfy));
dsdt = f_s+g_s*u;

%% autogen functions

params = [mt ml mw Iw l d lt r g];
refs = [xd hd];
gains = [kpx kph kdx kdh];

if ~exist('./2_link_gen')
    mkdir('./2_link_gen')
end
addpath('./2_link_gen')

matlabFunction(f_s, 'File', '2_link_gen/auto_f', 'Vars', {s,params}) ;
matlabFunction(g_s, 'File', '2_link_gen/auto_g', 'Vars', {s,params}) ;
matlabFunction(u, 'File', '2_link_gen/auto_u', 'Vars', {s,params,refs,gains}) ;
matlabFunction(dsdt, 'File', '2_link_gen/auto_dsdt', 'Vars', {s,params,refs,gains}) ;

matlabFunction(p1, 'File', '2_link_gen/auto_p1', 'Vars', {s,params}) ;
matlabFunction(p2, 'File', '2_link_gen/auto_p2', 'Vars', {s,params}) ;
matlabFunction(p3, 'File', '2_link_gen/auto_p3', 'Vars', {s,params}) ;
matlabFunction(p4, 'File', '2_link_gen/auto_p4', 'Vars', {s,params}) ;
matlabFunction(pe, 'File', '2_link_gen/auto_pe', 'Vars', {s,params}) ;

