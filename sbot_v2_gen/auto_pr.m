function pr = auto_pr(in1,in2)
%AUTO_PR
%    PR = AUTO_PR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Dec-2020 15:37:57

d = in2(:,5);
q1 = in1(4,:);
q2 = in1(5,:);
r = in2(:,6);
x = in1(1,:);
y = in1(2,:);
t2 = q1+q2;
t3 = 1.0./d;
t4 = r.*t2.*t3;
pr = [x+(d.*sin(t4))./2.0;y-(d.*cos(t4))./2.0;r];
