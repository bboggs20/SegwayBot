function draw_two_link_roller(s,params)

    r = params(8);

    p1 = auto_p1(s,params);
    p2 = auto_p2(s,params);
    p3 = auto_p3(s,params);
    p4 = auto_p4(s,params);

    pe = auto_pe(s,params);

    % Plot
    hold on
    rectangle('Position',[p1(1)-r 0 r*2 r*2],'Curvature',[1,1]); % wheel

    line([p1(1);p2(1)], [r;p2(2)], 'Color', 'k', 'LineWidth', 2); % leg
    line([p2(1);p3(1)], [p2(2);p3(2)], 'Color', 'k', 'LineWidth', 2); % torso

    plot(p2(1), p2(2), 'bo', 'MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g') % leg joint
    plot(p1(1), p1(2), 'bo', 'MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g') % wheel joint

    plot(pe(1), pe(2), 'bo', 'MarkerSize',3,'MarkerEdgeColor','g','MarkerFaceColor','b') % wheel encoder pos
    % plot(pCom(1), pCom(2), 'x', 'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5) % sys com pos
    % plot(pCom(1), pCom(2), 'bo', 'MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1.5) % sys com pos

    plot(p1(1), p1(2), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % wheel mass
    plot(p3(1), p3(2), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % torso mass
    plot(p4(1), p4(2), '*', 'MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m') % leg mass

end