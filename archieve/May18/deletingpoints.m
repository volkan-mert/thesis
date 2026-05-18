close all; clc

nichols(N_frd, 'r*')
ax = axes(gcf)
ax.CurrentPoint
ax = axes(gcf)
findobj(gcf,'Type','Line')
h = findobj(gcf,'Type','Line')
h.XData
h.XData(1)
h.YData
h.YData(1)
h
h(1)
delete(h(1))
h
h = findobj(gcf,'Type','Line')
h
h = findobj(gcf,'Type','Line')
nichols(N_frd, 'r*')

%%

h = findobj(gcf,'Type','Line')

xData = get(h,'XData')
xData = get(h,'YData')
xData = get(h,'XData')
yData = get(h,'YData')
xData(1) = []; yData(1) = [];
set(h,'XData',xData,'YData',yData);
h
drawnow

