
%% syms
clear

syms m mw Iw Iphi Ipsi g l r d x dx y dy thr dthr thl dthl dphi phi Kp Kd dthr_d dthl_d Tr Tl v real
sym psi real;

%% Lagrangian Dynamics

q = [phi; psi];
dq = [dphi; dpsi];

% psi = r/d*(thr-thl);
% dpsi = r/d*(dthr-dthl);
DR = [d/2 d/2; d/r -d/r];

TH = DR\[v; dpsi];
dthr = TH(1); dthl = TH(2);

%dx = -r/2*(dthr+dthl)*sin(psi);
%dy = r/2*(dthr+dthl)*cos(psi);

s = [x;y;v;q;dq];

Vc = v*[-sin(psi); cos(psi); 0] + [l*sin(phi)*sin(psi)*(dpsi-dphi); l*sin(phi)*cos(psi)*(dphi-dpsi); l*dphi*cos(phi)];
Wc = [-dphi*sin(phi)*sin(psi); dphi*sin(phi)*cos(psi); dphi*cos(phi) + dpsi];
I = [m*(l*sin(phi)*cos(psi))^2; m*(l*sin(phi)*sin(psi))^2; m*(l*cos(phi))^2];

KE = simplify(1/2*Iw*(dthr^2+dthl^2) + 1/2*mw*r^2*(dthr^2+dthl^2) + 1/2*Vc'*m*Vc + 1/2*I'*Wc.^2);
%KE = simplify((1/2*Iw + 1/2*mw*r^2)*((dthr-dphi)^2+(dthl-dphi)^2) + 1/2*Iphi*dphi^2 + 1/2*Ipsi*r^2/d^2*(dthr-dthl)^2 + 1/2*m*(dx^2+dy^2));
PE = m*g*l*cos(phi);

q_act = [thr; thl];

[D, C, G, B] = LagrangianDynamics(KE, PE, q, dq, q_act);
D = simplify(D);
C = simplify(C);
B = [-1 -1; 1 0; 0 1];
f_s = simplify([dq; D\(-C*dq-G)]);
g_s = simplify([zeros(3,2); D\B]);

%% functions for autogen

q = [phi; thr; thl];
dq = [dphi; dthr; dthl];

dpsi = dthr-dthl;
dpsid = dthr_d-dthl_d;
v = dthr+dthl;
vd = dthr_d+dthl_d;

s = [q;dq];

h = [dphi - (v-vd); dpsi-dpsid];
%h = [dphi; dthr-dthr_d; dthl-dthl_d];

Lfh = simplify(jacobian(h,s)*f_s);
Lgh = simplify(jacobian(h,s)*g_s);


%% autogen functions
if ~exist('./gen')
    mkdir('./gen')
end
addpath('./gen')

matlabFunction(f_s, 'File', 'gen/auto_f') ;
matlabFunction(g_s, 'File', 'gen/auto_g') ;
matlabFunction(Lgh, 'File', 'gen/auto_Lgh') ;
matlabFunction(Lfh, 'File', 'gen/auto_Lfh') ;


%% sim

clear m mw Iw Iphi Ipsi g l r d thr dthr thl dthl phi dphi psi Kp Kd dthr_d dthl_d Tr Tl z
close all;
addpath('./gen')

m = .150;
mw = .035;
r = .075;
Iw = mw*r^2;
g = 9.81;
l = .010;
d = .010;
Kp = 1; Kd = 0; Ki = 0;
dthr_d = 0;
dthl_d = 0;

tspan = [0 10];
s0 = [0;0;0;0;0;0];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t,s] = ode45(@sbot_dynamics_ode, tspan, s0, opts);

phi = s(:,1); thr = s(:,2); thl = s(:,3); 
dphi = s(:,4); dthr = s(:,5); dthl = s(:,6);

psi = (r/d).*(thr-thl);

dx = (-r/2).*(dthr+dthl).*sin(psi);
dy = (r/2).*(dthr+dthl).*cos(psi);
x = [dx(1)*t(1)]; y = [dy(1)*t(1)];
for i = 2:length(t)
    x = [x; dx(i)*(t(i)-t(i-1))];
    y = [y; dy(i)*(t(i)-t(i-1))];
end


% u = [];
% for i = 1:length(t)
%     Lfh = auto_Lfh(Iw,d,dphi(i),dthl(i),dthr(i),g,l,m,mw,phi(i),r,thl(i),thr(i));
%     Lgh = auto_Lgh(Iw,d,l,m,mw,phi(i),r,thl(i),thr(i));
%     
%     h = [dphi(i) - (dthr(i)+dthl(i)) + (dthr_d+dthl_d); dthr(i) - dthr_d - dthl(i) + dthl_d];
%     sh = [phi(i) - (thr(i)+thl(i)) + (dthr_d+dthl_d)*t(i); thr(i) - dthr_d*t(i) - thl(i) + dthl_d*t(i)];
%     v = (-Kp*h-Ki*sh)/(1+Kd);
%     
%     u = [u; (Lgh\(-Lfh+v))'];
% end

i = floor(length(x)/2);
figure(1); hold on; plot(x,y); plot_dir(x(i:i+1),y(i:i+1)); title("Trajectory"); xlabel("x [m]"); ylabel("y [m]"); xlim([-1 1]); ylim([-1 1]); hold off;
figure(2); hold on; plot(t,(mod(phi,2*pi)*180/pi)); title("\phi(t)"); xlabel("t [s]"); ylabel("\phi [deg]"); hold off;
figure(3); hold on; plot(t,dphi); title("d\phi/dt"); xlabel("t [s]"); ylabel("d\phi/dt"); hold off;
figure(4); hold on; plot(t,dthr,t,dthl); title("d\theta/dt"); xlabel("t [s]"); ylabel("d\theta/dt"); legend("\theta_r", "\theta_l"); hold off;
%figure(5); hold on; plot(t,u(:,1),t,u(:,2)); title("u(t)"); xlabel("t [s]"); ylabel("\tau [Nm]"); legend("\tau_r", "\tau_l"); hold off;

function dsdt = sbot_dynamics_ode(t,s)
    % constants (SI units)
    m = .150;
    mw = .035;
    r = .075;
    Iw = mw*r^2;
    g = 9.81;
    l = .010;
    d = .010;
    %Io = m*l^2;
    Kp = 1; Kd = 0; Ki = 0;
    dthr_d = 0;
    dthl_d = 0;
    %meu = 0.8;
    phi = s(1); thr = s(2); thl = s(3);
    dphi = s(4); dthr = s(5); dthl = s(6); 
    
    %psi = r/d*(thr-thl);
    %dpsi = r/d*(dthr-dthl);
    
    f_s = auto_f(Iw,d,dphi,dthl,dthr,g,l,m,mw,phi,r,thl,thr);
    g_s = auto_g(Iw,d,l,m,mw,phi,r,thl,thr);
    Lfh = auto_Lfh(Iw,d,dphi,dthl,dthr,g,l,m,mw,phi,r,thl,thr);
    Lgh = auto_Lgh(Iw,d,l,m,mw,phi,r,thl,thr);
    %h = [dphi; dthr-dthr_d; dthl-dthl_d];
    
    dpsi = dthr-dthl;
    dpsid = dthr_d-dthl_d;
    vc = dthr+dthl;
    vcd = dthr_d+dthl_d;
    
    h = [dphi - (vc-vcd); dpsi-dpsid];
    sh = [phi - (thr+thl-vcd*t); thr-thl-dpsid*t];
    v = -Kp*h;
    
    u = Lgh\(-Lfh+v);
    
    dsdt = f_s + g_s*u;
    fprintf("t = %f\n", t);
end







