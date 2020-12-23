function animate_two_link_roller(t,s,params)
    figure(1000)
    for i =1:length(t)-1 
        clf ;
        p1 = auto_p1(s(i, :)',params);
        x = p1(1);
        draw_two_link_roller(s(i, :)', params);
        line([-10, 20],[0;0],'Color', 'k', 'LineWidth', 2)
        xlim([x-.5 x+.5]);
        %ylim([-.25 3.5]);
        axis equal

        grid on ;
        drawnow ;
        pause(t(i+1)-t(i)) ;
    end
end