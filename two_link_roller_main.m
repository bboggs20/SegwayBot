
% system parameters
mt = .1;
ml = .1;
mw = .01;
Iw = .01;
l = .12;
d = .05;
lt = .2;
r = .035;
g = 9.81;

params = [mt ml mw Iw l d lt r g];

% reference output
xd = 5;
hd = .3;

refs = [xd hd];

% controller gains
kpx = 100;
kph = 100;
kdx = 5;
kdh = 5;

gains = [kpx kph kdx kdh];

% simulation
tspan = [0 10];
s0 = [pi/6;pi;0;0;0;0];
[t,s] = ode45(@(t,s)odefun(t,s,params,refs,gains), tspan, s0);

animate_two_link_roller(t,s,params);

% odefun
function sdot = odefun(t,s,params,refs,gains)
    f_s = auto_f(s,params);
    g_s = auto_g(s,params);
    u = auto_u(s,params,refs,gains);
    
    sdot = f_s+g_s*u;
end


