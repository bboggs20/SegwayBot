function draw_sbot_v2(s,params)
    addpath('./sbot_v2_gen')
    r = params(6);

    pb = auto_pb(s,params);
    %pr = auto_pr(s,params);
    %pl = auto_p3(s,params);
    pt = auto_pt(s,params);

    %pe = auto_pe(s,params);

    % Plot
    hold on
    rectangle('Position',[pb(1)-r 0 r*2 r*2],'Curvature',[1,1]); % wheel

    line([pb(1);pt(1)], [r;pt(3)], 'Color', 'k', 'LineWidth', 2); % leg
    %line([pr(1);pl(1)], [pr(2);pl(2)], 'Color', 'k', 'LineWidth', 2); % torso

    %plot(pr(1), pr(2), 'bo', 'MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g') % leg joint
    %plot(pb(1), pb(2), 'bo', 'MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g') % wheel joint

    %plot(pe(1), pe(2), 'bo', 'MarkerSize',3,'MarkerEdgeColor','g','MarkerFaceColor','b') % wheel encoder pos
    % plot(pCom(1), pCom(2), 'x', 'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5) % sys com pos
    % plot(pCom(1), pCom(2), 'bo', 'MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1.5) % sys com pos

    plot(pb(1), pb(3), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % wheel mass
    plot(pt(1), pt(3), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % torso mass
    %plot(pt(1), pt(2), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % leg mass

end