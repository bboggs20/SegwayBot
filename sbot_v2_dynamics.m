syms x y q1 q2 dq1 dq2 phi dphi m mw Iw l d r g real


% absolute coords
thr = q1 + phi;
thl = -q2 + phi;
dthr = dq1 + dphi;
dthl = -dq2 + dphi;

% yaw
psi = r/d*(thr-thl);
dpsi = r/d*(dthr-dthl);

% cartesian velocities
dx = r/2*(dthr+dthl)*cos(psi);
dy = r/2*(dthr+dthl)*sin(psi);

% relative coords
q = [phi;q1;q2]; 
dq = [dphi;dq1;dq2];
s = [x;y;q;dq];

% world positions
pb = [x; y; r];
pr = pb + d/2*[sin(psi); -cos(psi); 0];
pl = pb + d/2*[-sin(psi); cos(psi); 0];
pt = pb + l*[sin(phi)*cos(psi); sin(phi)*sin(psi); cos(phi)];

% world velocities
dpb = [dx; dy; 0];
dpr = dpb + d/2*dpsi*[cos(psi); sin(psi); 0];
dpl = dpb + d/2*dpsi*[-cos(psi); -sin(psi); 0];
dpt = dpb + l*dpsi*dphi*[cos(phi)*cos(psi) - sin(phi)*sin(psi); cos(phi)*sin(psi) + sin(phi)*cos(psi); -1/dpsi*sin(phi)];

% Lagrangian Dynamics
KE = .5*(dpr'*mw*dpr + dpl'*mw*dpl + dpt'*m*dpt + Iw*thr^2 + Iw*thl^2);
PE = g*(m*pt(3) + mw*(pr(3)+pl(3)));

qAct = [q1;q2]; % actuated angles

[D, C, G, B] = LagrangianDynamics(KE, PE, q, dq, qAct);

f_s = simplify([dx; dy; dq; D\(-C*dq-G)]);
g_s = simplify([zeros(2,size(qAct,1)); zeros(size(q,1),size(qAct,1)); D\B]);

%% Output linearization
syms phid psid kp kd real

h = [phi-phid; psi-psid];

Lfy = simplify(jacobian(h,s)*f_s);  % = dydt
Lgy = simplify(jacobian(h,s)*g_s);  % should be 0
LgLfy = simplify(jacobian(Lfy,s)*g_s);
Lf2y = simplify(jacobian(Lfy,s)*f_s);

u = simplify(pinv(LgLfy)*(-Lf2y -kp*h -kd*Lfy));
dsdt = f_s+g_s*u;

%% autogen
params = [m mw Iw l d r g];
gains = [kp kd];
refs = [phid psid];

if ~exist('./sbot_v2_gen')
    mkdir('./sbot_v2_gen')
end
addpath('./sbot_v2_gen')

matlabFunction(f_s, 'File', 'sbot_v2_gen/auto_f', 'Vars', {s,params}) ;
matlabFunction(g_s, 'File', 'sbot_v2_gen/auto_g', 'Vars', {s,params}) ;
matlabFunction(Lfy, 'File', 'sbot_v2_gen/auto_Lfy', 'Vars', {s,params}) ;
matlabFunction(LgLfy, 'File', 'sbot_v2_gen/auto_LgLfy', 'Vars', {s,params}) ;
matlabFunction(Lf2y, 'File', 'sbot_v2_gen/auto_Lf2y', 'Vars', {s,params}) ;
matlabFunction(psi, 'File', 'sbot_v2_gen/auto_psi', 'Vars', {s,params}) ;

matlabFunction(u, 'File', 'sbot_v2_gen/auto_u', 'Vars', {s,params,refs,gains}) ;
matlabFunction(dsdt, 'File', 'sbot_v2_gen/auto_dsdt', 'Vars', {s,params,refs,gains}) ;

matlabFunction(pb, 'File', 'sbot_v2_gen/auto_pb', 'Vars', {s,params}) ;
matlabFunction(pr, 'File', 'sbot_v2_gen/auto_pr', 'Vars', {s,params}) ;
matlabFunction(pl, 'File', 'sbot_v2_gen/auto_pl', 'Vars', {s,params}) ;
matlabFunction(pt, 'File', 'sbot_v2_gen/auto_pt', 'Vars', {s,params}) ;
