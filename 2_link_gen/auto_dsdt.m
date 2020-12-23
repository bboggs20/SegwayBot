function dsdt = auto_dsdt(in1,in2,in3,in4)
%AUTO_DSDT
%    DSDT = AUTO_DSDT(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Dec-2020 19:24:08

Iw = in2(:,4);
d = in2(:,6);
dq1 = in1(4,:);
dq2 = in1(5,:);
dq3 = in1(6,:);
l = in2(:,5);
lt = in2(:,7);
ml = in2(:,2);
mt = in2(:,1);
mw = in2(:,3);
q1 = in1(1,:);
q2 = in1(2,:);
r = in2(:,8);
t2 = cos(q1);
t3 = cos(q2);
t4 = q1+q2;
t5 = d.^2;
t6 = l.^2;
t7 = lt.^2;
t8 = ml.^2;
t9 = mt.^2;
t10 = q1.*2.0;
t11 = q2.*2.0;
t12 = r.^2;
t13 = Iw.*lt.*2.0;
t17 = 1.0./lt;
t19 = 1.0./mt;
t20 = -q2;
t14 = cos(t10);
t15 = cos(t11);
t16 = cos(t4);
t18 = 1.0./t7;
t21 = q2+t4;
t22 = q1+t4;
t25 = lt.*mt.*t12;
t26 = q1+t20;
t27 = Iw.*mt.*t6;
t28 = Iw.*l.*t3.*2.0;
t29 = t4.*2.0;
t31 = lt.*ml.*t12.*2.0;
t32 = lt.*mw.*t12.*2.0;
t33 = Iw.*ml.*t5.*2.0;
t35 = d.*l.*ml.*r.*t2;
t38 = d.*ml.*t3.*t12;
t39 = d.*l.*ml.*mt.*t12;
t40 = mt.*r.*t2.*t6;
t41 = l.*mt.*t3.*t12;
t43 = ml.*mt.*t5.*t12;
t44 = ml.*r.*t2.*t5.*2.0;
t45 = ml.*mt.*t6.*t12;
t46 = mt.*mw.*t6.*t12;
t48 = l.*ml.*t3.*t12.*2.0;
t49 = l.*mw.*t3.*t12.*2.0;
t53 = t5.*t8.*t12;
t54 = t6.*t9.*t12;
t56 = ml.*mw.*t5.*t12.*2.0;
t23 = cos(t21);
t24 = cos(t22);
t30 = cos(t26);
t34 = t27.*2.0;
t36 = -t28;
t37 = cos(t29);
t42 = l.*lt.*mt.*r.*t16;
t47 = t39.*2.0;
t50 = t15.*t27;
t51 = t14.*t25;
t52 = d.*lt.*ml.*r.*t16.*2.0;
t55 = t43.*2.0;
t57 = -t44;
t58 = t45.*2.0;
t59 = t46.*2.0;
t60 = -t39;
t62 = -t48;
t63 = -t40;
t64 = -t41;
t65 = -t49;
t72 = t14.*t39;
t73 = t15.*t39;
t77 = t14.*t43;
t78 = t15.*t45;
t79 = t15.*t46;
t61 = -t47;
t66 = d.*l.*ml.*r.*t23;
t67 = -t52;
t68 = -t42;
t69 = d.*ml.*t12.*t24;
t70 = l.*mt.*t12.*t24;
t71 = mt.*r.*t6.*t23;
t74 = -t50;
t75 = l.*lt.*mt.*r.*t30;
t76 = -t51;
t80 = -t77;
t81 = -t78;
t82 = -t79;
t83 = t37.*t39;
t84 = t37.*t47;
t85 = t37.*t53;
t86 = t37.*t54;
t87 = t37.*t60;
t88 = t83.*-2.0;
t89 = -t85;
t90 = -t86;
t91 = t13+t25+t31+t32+t35+t36+t38+t57+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t75+t76;
t92 = t27+t33+t43+t45+t46+t53+t56+t60+t72+t73+t74+t80+t81+t82+t87+t89;
t93 = 1.0./t92;
t94 = t17.*t91.*t93;
t95 = -t94;
dsdt = [dq1;dq2;dq3;(-t17.*t93.*(t28-t35-t38+t40+t41+t44+t48+t49-t66-t69-t70-t71)).*Inf+(-t18.*t19.*t93.*(t33+t34+t53+t54+t55+t56+t58+t59+t61+t88+t89+t90+d.*ml.*t3.*t25+d.*ml.*t24.*t25-l.*ml.*t3.*t25.*2.0-l.*mw.*t3.*t25.*2.0-Iw.*l.*lt.*mt.*t3.*2.0-l.*lt.*t3.*t9.*t12+l.*lt.*t9.*t12.*t24)).*Inf;t95.*Inf+(t18.*t19.*t93.*(t33+t34+t53+t54+t55+t56+t58+t59+t61+t88+t89+t90+Iw.*mt.*t7.*2.0+t7.*t9.*t12+d.*ml.*t3.*t25.*2.0+d.*ml.*t24.*t25.*2.0-l.*ml.*t3.*t25.*4.0-l.*mw.*t3.*t25.*4.0+ml.*mt.*t7.*t12.*2.0+mt.*mw.*t7.*t12.*2.0-t7.*t9.*t12.*t14-Iw.*l.*lt.*mt.*t3.*4.0-l.*lt.*t3.*t9.*t12.*2.0+l.*lt.*t9.*t12.*t24.*2.0)).*Inf;t95.*Inf+(t93.*(Iw.*2.0+ml.*t5.*2.0+ml.*t12.*2.0+mt.*t6+mt.*t12+mw.*t12.*2.0-mt.*t6.*t15-mt.*t12.*t14-d.*ml.*r.*t16.*4.0-l.*mt.*r.*t16.*2.0+l.*mt.*r.*t30.*2.0)).*Inf];
