
animate = 0; % 1 to animate

% system parameters
m = .1;
mw = .01;
Iw = .01;
l = .12;
d = .05;
r = .035;
g = 9.81;

params = [m mw Iw l d r g];

% reference outputs
phid = deg2rad(3);
psid = deg2rad(-1);
refs = [phid psid];

% controller gains
kp = 500; kd = 50;
gains = [kp kd];

% simulation
f = 50; %hz

tspan = (0:1/f:10);
s0 = [0;0;0;0;0;-3;0;0];
[t,s] = ode45(@(t,s)odefun(t,s,params,refs,gains), tspan, s0);

x = s(:,1); y = s(:,2); phi = s(:,3);

psi = zeros(size(t));
u = zeros(length(t),2);
for i=1:length(t)
    psi(i) = rad2deg(auto_psi(s(i,:)',params));
    
    kp = gains(1); kd = gains(2);
    Lfy = auto_Lfy(s(i,:)',params);
    LgLfy = auto_LgLfy(s(i,:)',params);
    Lf2y = auto_Lf2y(s(i,:)',params);
    y = s(i,3) - phid;
    u(i,:) = pinv(LgLfy)*(-Lf2y - kp*y - kd*Lfy);
end

figure(1);
plot(t,psi);

figure(2);
plot(x,y);

figure(3);
plot(t,rad2deg(phi));

figure(4); 
plot(t,u(:,1),t,-u(:,2));
legend('u_r','-u_l');

if animate == 1
    animate_sbot_v2(t,s,params);
end

% odefun
function sdot = odefun(t,s,params,refs,gains)
%     kp = gains(1); kd = gains(2);
% 
%     f_s = auto_f(s,params);
%     g_s = auto_g(s,params);
%     Lfy = auto_Lfy(s,params);
%     LgLfy = auto_LgLfy(s,params);
%     Lf2y = auto_Lf2y(s,params);
%     y = s(3) - phid;
%     
%     u = LgLfy\(-Lf2y - kp*y - kd*Lfy);
    
    sdot = auto_dsdt(s,params,refs,gains);
end


