function pb = auto_pb(in1,in2)
%AUTO_PB
%    PB = AUTO_PB(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-Dec-2020 15:37:57

r = in2(:,6);
x = in1(1,:);
y = in1(2,:);
pb = [x;y;r];
