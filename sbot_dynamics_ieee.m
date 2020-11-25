
%% symbolics

syms x y psi phi dphi v dpsi m mw r l d g Iwa Ix Ipsi Iphi Iwd kqv kq kav ka psi_d phi_d phi_m kr kv k_ v_d z7 real

Dphi = (m*r*l*cos(phi))^2 - r^2*(l^2*(m^2+2*m*mw)+2*Iphi*mw+Iphi*m) - 2*m*l^2*Iwa - 2*Iphi*Iwa;
Gphi = r^2*cos(phi)^2*(Ipsi-Ix-2*m*l^2*Iwa) + r^2*(m*l^2+Ix+2*Iwd+d^2/2*mw) + d^2/2*Iwa;

g_x = [zeros(4,2);
    ((m*r^2 + 2*mw*r^2 + 2*Iwa + m*r*l*cos(phi))/Dphi)*ones(1,2);
    (-(r*(m*r*l*cos(phi) + Iphi + m*l^2))/Dphi)*ones(1,2);
    r*d/(2*Gphi) -r*d/(2*Gphi)];

H_ = (3/2*m*r^2 + Iwa)*(Ipsi-Ix) - m*l^2*(Iwa + mw*r^2);
Kphi = (sin(phi)*(m*r^2*l*(Ix-Ipsi) - 4*Iphi*m*r^2*l - 3*m^2*r^2*l^3) + sin(3*phi)*(m*r^2*l*(Ix-Ipsi) + m^2*r^2*l^3))/(4*Dphi);

f1 = [v*cos(psi);
      v*sin(psi);
      dpsi;
      dphi];

f2 = [H_*dpsi^2*sin(2*phi)/Dphi + (m*r*l*dphi)^2*sin(2*phi)/(2*Dphi) - g*sin(phi)*(m^2*r^2*l + 2*Iwa*m*l + 2*m*r^2*l)/Dphi;
      Kphi*dpsi^2 + (m*r*l)^2*g*sin(2*phi)/(2*Dphi) - dphi^2*sin(phi)*(Iphi*m*r^2*l + m^2*r^2*l^3)/Dphi;
      (-r^2*dpsi*dphi*sin(2*phi)*(m*l^2 + Ix + Ipsi) - m*r^2*l*v*dpsi*sin(phi))/Gphi];

f = [f1; f2];

f21phi = - g*sin(phi)*(m^2*r^2*l + 2*Iwa*m*l + 2*m*r^2*l)/Dphi;
f22phi = (m*r*l)^2*g*sin(2*phi)/(2*Dphi);

dvss = simplify(1/g_x(5,1)*(g_x(5,1)*f22phi - f21phi*g_x(6,1)));
dvssphi = subs(dvss, phi, phi_d);

dphi_d = simplify((-kr*dvss-kv*(phi_m^2-phi_d^2)^2*(v-v_d)*dvssphi/phi_d)*exp(-k_*abs(dphi)));

w1 = -kqv*dpsi-kq*(psi-psi_d);
w2 = -kav*dphi-ka*(phi-phi_d);

%z7 = -dphi*g_x(6,1)+v*g_x(5,1);

z = [psi;dpsi;phi;dphi;x;y;z7;phi_d];



zz7 = -dphi*g_x(6,1)+v*g_x(5,1);
vz = solve(z7 == zz7,v);

dzz = [dpsi;w1;dphi;w2;
      ((z7+g_x(6,1)*dphi)/g_x(5,1))*cos(psi);
      ((z7+g_x(6,1)*dphi)/g_x(5,1))*sin(psi);
      dphi*(-dphi*jacobian(g_x(6,1),phi) + jacobian(g_x(5,1),phi)*(z7+g_x(6,1)*dphi)/g_x(5,1)) + g_x(5,1)*(f2(2)-f2(1)*g_x(6,1)/g_x(5,1));
      dphi_d];

dz = subs(dzz, v, vz);
%% constants


% phi_m: phi (pitch) angle bound (+-)
% kr, kv, k_: outer loop gains (phi_d gains)
    % k_ controls virtual time delta on new phi_d outputs
% kqv, kq: psi (yaw) D and P gains respectively
% kav, ka: phi (pitch) D and P gains, respectively


params = [m, mw, r, l, d, g, Iwa, Ix, Ipsi, Iphi, Iwd, phi_m];
gains = [kqv, kq, kav, ka, kr, kv, k_];

% SI units
pvals = [0.1, 0.01181, 0.0305, 0.07, 0.094, 9.81, 0.00000361, 0.0000001, 0.0000001, 0.000245, 0.00000712, .707];
kvals = [10,100,110,3000,55,5,2000];


%% autogen functions




dz_ = simplify(subs(dz, [params,gains], [pvals,kvals]));

if ~exist('./gen')
    mkdir('./gen')
end
addpath('./gen')

matlabFunction(dz_, 'File', 'gen/auto_dz', 'Vars', {z, psi_d, v_d}) ;
matlabFunction(vz, 'File', 'gen/auto_v', 'Vars', {params, gains, z7, phi, dphi}) ;


%% sim


tspan = [0 25];
s0 = [0;0;pi/288;0;0;0;-7000;0.0001];
psi_d = 0; v_d = 0;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t,z] = ode45(@(t,z)sbot_ode(t,z,psi_d,v_d), tspan, s0, opts);

phi_rec = z(:,3).*180./pi;


v_rec = zeros(size(t));
for i=1:length(t)
    v_rec(i) = auto_v(pvals, kvals, z(i,7), z(i,3), z(i,4));
    fprintf("t = %f s\n",t(i));
end

phi_d_rec = z(:,8).*180./pi;

x_rec = z(:,5); y_rec = z(:,6);

figure(1); plot(t,phi_rec,t,phi_d_rec); xlabel("t [s]"); ylabel("\phi [deg]"); legend("\phi","\phi_d");
figure(2); plot(t,v_rec); xlabel("t [s]"); ylabel("v [m/s]");
figure(3); plot(x_rec,y_rec); title("Trajectory"); xlabel("x [m]"); ylabel("y [m]");




%% odefun
function dzdt = sbot_ode(t,z,psi_d,v_d)
    dzdt = auto_dz(z,psi_d,v_d);
    
    fprintf("t = %f s\n",t);
end




















