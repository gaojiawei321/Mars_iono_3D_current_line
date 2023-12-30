function [outputArg1,outputArg2] = plot_3d_mag_global_quiv(qx,Xgrid,Ygrid,Zgrid,LString,qx1,qy,qz)
%PLOT_3D_MAG_GLOBAL Summary of this function goes here
%   Detailed explanation goes here



%yslice = -1.5:0.75:1.5;
xslice = -1.5:0.5:1.5;
%xslice = -0:0.25:1.5;

yslice = [];
zslice = [];

%yslice = [-1.0000        0   1.0000   ];
%zslice = [-1.0000        0   1.0000  ];

ys=3;
zs=3;

xsj=3;

v1=-9;
v2=12;

Cjx=zeros(length(Xgrid),length(Ygrid),length(Zgrid));
Cjx(1:end-1,1:end-1,1:end-1)=qx;
Cjy=zeros(length(Xgrid),length(Ygrid),length(Zgrid));
Cjy(1:end-1,1:end-1,1:end-1)=qy;
Cjz=zeros(length(Xgrid),length(Ygrid),length(Zgrid));
Cjz(1:end-1,1:end-1,1:end-1)=qz;

[Xd,Yd,Zd]=meshgrid(Xgrid,Ygrid,Zgrid);

%slice(X,Y,Z,curlx, xslice,[],[]); hold on;

slice(Xd,Yd,Zd, Cjx, xslice,[],[]); hold on;


%[startX,startY,startZ] = meshgrid(-1.5,-2:0.2:2,-2:0.2:2);
%[startX,startY,startZ] = meshgrid([-1,0],-2:0.2:2,-2:0.2:2);


%verts = stream3(X,Y,Z,curlx,curly,curlz,startX,startY,startZ);
%streamline(verts); hold on;

% quiver in different slice to be test...


xjg=(Xgrid(1:end-1)+Xgrid(2:end))/2;
yjg=(Xgrid(1:end-1)+Xgrid(2:end))/2;
zjg=(Xgrid(1:end-1)+Xgrid(2:end))/2;


%[qux,quy,quz] = meshgrid(Xgrid,Xgrid,Xgrid);
[qux,quy,quz] = meshgrid(xjg,yjg,zjg);

%[Xq,Yq,Zq] = meshgrid(xslice,yjg,zjg);
[Xq,Yq,Zq] = meshgrid(xslice,-2.8:0.4:2.8,-2.8:0.4:2.8);

%yjg=-2.8:0.4:2.8;
%zjg=-2.8:0.4:2.8;


qx1=interp3(qux,quy,quz,qx,Xq,Yq,Zq);
qy1=interp3(qux,quy,quz,qy,Xq,Yq,Zq);
qz1=interp3(qux,quy,quz,qz,Xq,Yq,Zq);


%quiver3(Xq,Yq,Zq,qx1*0,qy1,qz1,5,'k'); hold on;
hq=quiver3(Xq,Yq,Zq,qx1*0,qy1,qz1,7,'k','LineWidth',1,'MaxHeadSize',0.02,'AutoScale','on'); hold on;
%hq=quiver3(Xq,Yq,Zq,qx1*0,qy1,qz1,6,'k','LineWidth',1,'MaxHeadSize',0.02); hold on;



%text(-2.26,0.116951717138591,1.335514753144047,'\rightarrow','Rotation',90,'FontSize',25); hold on;

%text(-2.193,0.117202713050204,0.995194374340791,'10 nT','Rotation',90,'FontSize',15); hold on;

%text(-2.26,0.11,2.23,'\rightarrow','Rotation',90,'FontSize',25); hold on;
%text(-2.19,0.11,2.18,'10 nT','Rotation',90,'FontSize',12); hold on;


%'MarkerSize',4

%for kj=4%1:length(xslice)
%     xkj=xslice(kj);
%quiver3(qux(30,:,:),quy,quz,qx(30,:,:),qy,qz,'k'); hold on;
%end

%quiver3(X,Y,Z,curlx,curly,curlz,'k');


for kj=1:length(xslice)
    xkj=xslice(kj);
    plot3([xkj,xkj,xkj,xkj,xkj],[-xsj,xsj,xsj,-xsj,-xsj],[-xsj,-xsj,xsj,xsj,-xsj],'k','LineWidth',1) ;hold on;
    
    %bow shock model Nˇemec
    a=4.219;c=1.464;b=-0.063;gamma=0.205;delta=0.018;
    Psw=0.637; Nsw=2.34;Vsw=392;B=4.45; F=1.087;
    x0bs=c*Psw^b*F^gamma*B^0.018;
    p_bs=(a*(x0bs-xkj)).^0.5;
    
    
    %bow shock model Nˇemec
    a=1.567;c=1.187;b=-0.065;gamma=0.094;delta=0.038;
    Psw=0.637; Nsw=2.34;Vsw=392;B=4.45; F=1.087;
    x0mpb=c*Psw^b*F^gamma*B^0.018;
    p_mpb=(a*(x0mpb-xkj)).^0.5;
    
    theta=0:0.01:2*pi;
    y= 0+cos(theta);
    z= 0+sin(theta); % plot model bs and mpb
    plot3(xkj*ones(length(y),1),y*p_bs,z*p_bs,'LineStyle','--','color','r','linewidth',2); hold on;
    plot3(xkj*ones(length(y),1),y*p_mpb,z*p_mpb,'LineStyle','--','color','m','linewidth',2); hold on;
    
    % plot Mars body
    
    Rm1=(1-(xkj)^2)^0.5;
     plot3(xkj*ones(length(y),1),y*Rm1,z*Rm1,'k','linewidth',2); hold on;
    
end




xl1=xlabel('X_{MSE} (R_m)','FontSize',15);

%xl1.Position=[xl1.Position(1)+0.5 xl1.Position(2) xl1.Position(3)];

yl1=ylabel('Y_{MSE} (R_m)','FontSize',15);


zlabel('Z_{MSE} (R_m)','FontSize',15)

%zlabel('Z_{MSE} (R_m)','FontWeight','bold')

caxis([-10,10]);
h1=colorbar;
h1.Label.String = [ LString];
h1.FontSize=15;
h1.TickDirection='out';
h1.FontSize=15;
h1.TickDirection='out';
 h1.TickLength=0.05;
h1.LineWidth=1

%h1.FontAngle='italic';
h1.FontName='times';
 
%h1.Label.Position = [ 3.490624904632568,-3.276826739311218,0];
% arrow and text same line


shading flat
%title('Bcrust<10 nT  J_x')
% title('Bcrust<10 nT  B_{\theta} ')
% title('Bcrust<10 nT  B_{\rho} std')


set(gca,'YDir','reverse');
set(gca,'XDir','reverse');

%view(76,13)
view(v1,v2)

xlim([-1.5 1.5])

ylim([-ys ys])
zlim([-zs zs])

%colormap(mycolormap7)
set(gca,'linewidth',1.5)

set(gca,'Fontsize',12);

xl1.Position=[ 0.2   3.   -4.];
yl1.Position=[1.60 0.2 -3.5];

a=1

end

