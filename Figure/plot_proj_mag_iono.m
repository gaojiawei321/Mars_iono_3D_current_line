function [outputArg1,outputArg2] = plot_proj_mag_iono(qx,Ygrid,Zgrid,qr,qt,qp)
%PLOT_PROJ_MAG_IONO Summary of this function goes here
%   Detailed explanation goes here

ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
    axesm mollweid ;
    axis off
    dc=0;
framem; %gridm
gridm


xx1 =squeeze(qx);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

k = (1/4)*ones(2);

xx2=conv2(xx1,k,'same');
xx2(end,:)=xx1(end,:);
xx2(:,end)=xx1(:,end);

qt=conv2(qt,k,'same');
qp=conv2(qp,k,'same');

%qt=mapminmax(qt);
%qp=mapminmax(qp);


% add addtional line
xx1j=[xx2,zeros(sz1(2),1)];
xx1j=[xx1j;zeros(1,sz1(3)+1)];


pcolorm(Ygrid/pi*180,Zgrid/pi*180,xx1j)

[gr1,gr2]=meshgrid(ygj/pi*180,zgj/pi*180);
q1=quiverm(gr1,gr2,-qt',qp'); hold on;

q1(1).Color=[0.1 0.1 0.1 ];
q1(2).Color=[0.1 0.1 0.1 ];

q1(1).LineWidth=2;
q1(2).LineWidth=2;
%contourfm(Ygrid/pi*180,Zgrid/pi*180,xx1j,-20:5:20)


ter1=90*ones(1,100);
ter2=-90*ones(1,100);
ter3=0*ones(1,100);

tlat=linspace(-90,90,100);

plotm(tlat,ter1,'Linestyle','--','Color','k','Linewidth',1.5)
plotm(tlat,ter2,'Linestyle','--','Color','k','Linewidth',1.5)
plotm(tlat,ter3,'Linestyle','--','Color','k','Linewidth',1.5)


%% text
%textm(-90,-180,'-','FontSize',10,'fontweight', 'normal','HorizontalAlignment','left','VerticalAlignment','middle')
%textm(90,-180,'-','FontSize',10,'fontweight', 'normal','HorizontalAlignment','right','VerticalAlignment','middle')

%textm(-45,180,'-','FontSize',10,'fontweight', 'normal','HorizontalAlignment','left')
%textm(0,180,'��','FontSize',10,'fontweight', 'normal','HorizontalAlignment','left','VerticalAlignment','middle')
%textm(45,180,'��','FontSize',10,'fontweight', 'normal','HorizontalAlignment','left')


text(2.94,0,'0\circ','FontSize',10,'fontweight', 'normal')
text(2.54,0.84,'45\circ  N','FontSize',10,'fontweight', 'normal')
text(2.54,-0.84,'45\circ S','FontSize',10,'fontweight', 'normal')



text(-3.2,0,'0\circ','FontSize',10,'fontweight', 'normal')
text(-3.00,0.84,'45\circ N','FontSize',10,'fontweight', 'normal')
text(-3.00,-0.84,'45\circ S','FontSize',10,'fontweight', 'normal')
%text(-3.5,0.84,'45\circ N','FontSize',10,'fontweight', 'normal')
%text(-3.5,-0.84,'45\circ S','FontSize',10,'fontweight', 'normal')


%S9, S10 yasuo 
%text(-3.39,0,'0\circ','FontSize',10,'fontweight', 'normal')
%text(-3.28,0.84,'45\circ','FontSize',10,'fontweight', 'normal')
%text(-3.28,-0.84,'45\circ','FontSize',10,'fontweight', 'normal')



%text(0,1.54,'90\circ N','FontSize',10,'fontweight', 'normal')
%text(0,-1.54,'90\circ S','FontSize',10,'fontweight', 'normal')

%text(1.5,0,'90\circ E','FontSize',10,'fontweight', 'normal','VerticalAlignment','top')
%text(-1.3,0,'90\circ W','FontSize',10,'fontweight', 'normal','VerticalAlignment','top')

%text(1.4,0,'18','FontSize',10,'fontweight', 'normal','VerticalAlignment','top')
%text(-1.4,0,'6','FontSize',10,'fontweight', 'normal','VerticalAlignment','top')
%text(0,0,'12','FontSize',10,'fontweight', 'normal','VerticalAlignment','top')


%%


%q1=quiverm(0,45,10,10,'r'); hold on;



%shading interp

%colormap(mycolormap8)


caxis([-30,30]);
%h1=colorbar;
%h1.Label.String = '(nA/m^2)';

 



end

