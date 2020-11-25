T_s = .0785; % Nm
w_nl = 140; % rpm
w = (0:.01:w_nl);

m = -T_s/w_nl;

T = m.*w + T_s;

figure; hold on; grid on;

plot(w,T);
title("Torque-Speed Curve");
xlabel("\omega [rpm]");
ylabel("\tau [Nm]");