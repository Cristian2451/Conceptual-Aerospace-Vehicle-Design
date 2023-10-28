clear
clc
s = settings;

x(16) = 0; % initialise
x(1) = 0;

ta = 200; % turnaround
x(2) = 100;

x(3) = x(2) + 79.61; % leg 1 climb
x(4) = x(3) + 200.78; % leg 1 cruise
x(5) = x(4) + 79.61; % leg 1 descent

x(6) = x(5) + ta;

x(7) = x(6) + 94.09; % leg 2 climb
x(8) = x(7) + 521.82; % leg 2 cruise
x(9) = x(8) + 94.09; % leg 2 descent

x(10) = x(9) + ta;

x(11) = x(10) + 102; % leg 3 climb
x(12) = x(11) + 219.2; % leg 3 cruise
x(13) = x(12) + 98.8; % leg 3 descent

x(14) = x(13) + 102; % diversion climb
x(15) = x(14) + 178.79; % diversion cruise
x(16) = x(15) + 89.21; % diversion descent

x(17) = x(16) + 300; % loiter
x(18) = x(17) + 15.99; % final descent


h1 = 25000; % cruise 1 altitude
h2 = 29500; % cruise 2 altitude
h3 = 32000; % cruise 3 altitude
hmiss = 1000; % missed approach min altitude
hdiv = 33000; % diversion altitude
hloi = 5000; % loiter altitude

y = [0,0,h1,h1,0,0,h2,h2,0,0,h3,h3,hmiss,hdiv,hdiv,hloi,hloi,0];

% plot up to loiter
plot(x(1:16),y(1:16),'color',[0.4940 0.1840 0.5560],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.4940 0.1840 0.5560])
hold on
% plot loiter
plot(x(16:17),y(16:17),'--','color',[0.4940 0.1840 0.5560],'LineWidth',2)
% plot after loiter
plot(x(17:18),y(17:18),'color',[0.4940 0.1840 0.5560],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.4940 0.1840 0.5560])

set(gca,'xtick',[])
ylabel('Altitude / ft','FontName','Verdana','FontWeight','bold','FontSize',13)
xlabel({'';'';'Approximate Relative Mission Duration'},'FontName','Verdana','FontWeight','bold','FontSize',13)
ylim([0,37000])
set(gca,'YTickLabel',get(gca,'YTick'),'FontSize',12)
grid on

% leg labels

text((x(1)+x(2))/2,0,{'\uparrow','Warmup & Taxi'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(2)+x(3))/2,h1/2,'\leftarrow Climb 1','FontName','Verdana','FontSize',12)
text((x(3)+x(4))/2,h1,{'Cruise 1','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(4)+x(5))/2,h1/2,'\leftarrow Descent 1','FontName','Verdana','FontSize',12)
text((x(5)+x(6))/2,0,{'\uparrow','Christchurch Turnaround'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(6)+x(7))/2,h2/2,'\leftarrow Climb 2','FontName','Verdana','FontSize',12)
text((x(7)+x(8))/2,h2,{'Cruise 2','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(8)+x(9))/2,h2/2,'\leftarrow Descent 2','FontName','Verdana','FontSize',12)
text((x(9)+x(10))/2,0,{'\uparrow','Hamilton Turnaround'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(10)+x(11))/2,h3/2,'\leftarrow Climb 3','FontName','Verdana','FontSize',12)
text((x(11)+x(12))/2,h3,{'Cruise 3','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((2*x(12)+x(13))/3,(2*h3-hmiss)/3,'Descent 3 \rightarrow','HorizontalAlignment','right','FontName','Verdana','FontSize',12)
text((9*x(13)+x(14))/10,(hdiv+9*hmiss)/10,'\leftarrow Diversion Climb','FontName','Verdana','FontSize',12)
text((x(14)+x(15))/2,hdiv,{'Diversion Cruise','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(15)+x(16))/2,(hdiv+hloi)/2,'\leftarrow Diversion Descent','FontName','Verdana','FontSize',12)
text((x(16)+x(17))/2,hloi,{'Loiter','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana','FontSize',12)
text((x(17)+x(18))/2,hloi/2,'\leftarrow Final Descent','FontName','Verdana','FontSize',12)

hold off



% % numbers instead
% 
% figure
% 
% plot(x,y,'color',[0.4940 0.1840 0.5560],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.4940 0.1840 0.5560])
% set(gca,'xtick',[])
% ylabel('Altitude / ft','FontName','Verdana','FontWeight','bold')
% xlabel({'';'';'Approximate Relative Mission Duration'},'FontName','Verdana','FontWeight','bold')
% ylim([0,40000])
% set(gca,'YTickLabel',get(gca,'YTick'))
% 
% % leg labels
% 
% text((x(1)+x(2))/2,0,{'\uparrow','1'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana')
% text((x(2)+x(3))/2,h1/2,'\leftarrow 2','FontName','Verdana')
% text((x(3)+x(4))/2,h1,{'3','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana')
% text((x(4)+x(5))/2,h1/2,'\leftarrow 4','FontName','Verdana')
% text((x(5)+x(6))/2,0,{'\uparrow','5'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana')
% text((x(6)+x(7))/2,h2/2,'\leftarrow 6','FontName','Verdana')
% text((x(7)+x(8))/2,h2,{'7','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana')
% text((x(8)+x(9))/2,h2/2,'\leftarrow 8','FontName','Verdana')
% text((x(9)+x(10))/2,0,{'\uparrow','9'},'VerticalAlignment','top','HorizontalAlignment','center','FontName','Verdana')
% text((x(10)+x(11))/2,h3/2,'\leftarrow 10','FontName','Verdana')
% text((x(11)+x(12))/2,h3,{'11','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana')
% text((2*x(12)+x(13))/3,(2*h3-hmiss)/3,'12 \rightarrow','HorizontalAlignment','right','FontName','Verdana')
% text((2*x(13)+x(14))/3,(hdiv+2*hmiss)/3,'\leftarrow 13','FontName','Verdana')
% text((x(14)+x(15))/2,hdiv,{'14','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana')
% text((x(15)+x(16))/2,(hdiv+hloi)/2,'\leftarrow 15','FontName','Verdana')
% text((x(16)+x(17))/2,hloi,{'16','\downarrow'},'VerticalAlignment','bottom','HorizontalAlignment','center','FontName','Verdana')
% text((x(17)+x(18))/2,hloi/2,'\leftarrow 17','FontName','Verdana')

