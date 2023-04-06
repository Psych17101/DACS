set(0,'DefaultFigureVisible','off');

firstfig = hgload('n=1_nuxy_inpl.fig');
secondfig = hgload('n=2_nuxy_inpl.fig');
thirdfig = hgload('n=3_nuxy_inpl.fig');

figure(20)
h(1) = subplot(1,1,1);

firstplot = copyobj(allchild(get(firstfig,'CurrentAxes')),h(1));
set(firstplot,'LineStyle', '-', 'Color', 'r');
hold on;
secondplot = copyobj(allchild(get(secondfig,'CurrentAxes')),h(1));
set(secondplot,'LineStyle', '--', 'Color', 'b');
thirdplot = copyobj(allchild(get(thirdfig,'CurrentAxes')),h(1));
set(thirdplot,'LineStyle', ':', 'Color', 'k');
hold off;

xlabel('θ (degrees)');
ylabel('ν_x_y');
set(gca,'FontSize',14)
grid on;

l(1) = legend([firstplot, secondplot, thirdplot],'n=1','n=2',...
    'n=3');


set(20,'Visible','on');