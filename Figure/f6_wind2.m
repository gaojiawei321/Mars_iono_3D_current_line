
clear all;close all;clc

load ionMSOsph_Brtp_L10_Bm120_2023

%load ionMSOsph_Brtp_L10_Bm120_2023_dLS180


%load ionMSOsph_Brtp_L10_Bm120_Ls0_90
%load ionMSOsph_Brtp_L10_Bm120_Ls90_180
%load ionMSOsph_Brtp_L10_Bm120_Ls270_360
%load ionMSOsph_Brtp_L10_Bm120_Ls270_360


load mycolormap9



%%

Rm=3393.5;
[xge,yge,zge]=size(quanBr);

%xge=10;
%yge=10;
%zge=20;

%xs=100+Rm;
%xe=600+Rm;

%ys=-pi/2;
%ye=pi/2;

%zs=-pi;
%ze=pi;

%xj=(xe-xs)/xge;
%yj=(ye-ys)/yge;
%zj=(ze-zs)/zge;

%Xgrid=linspace(xs,xe,xge+1);
%Ygrid=linspace(ys,ye,yge+1);
%Zgrid=linspace(zs,ze,zge+1);

xgj=Xgrid(1:end-1)+((Xgrid(2)-Xgrid(1)))/2;  % r grid
xgj=xgj*1000; % km --> m
ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

ge=zeros(xge,yge,zge);



%%     using sph curl

Br=quanBr.*10^-9;
Bt=quanBt.*10^-9;
Bp=quanBp.*10^-9;


for i=1:xge
    for j=1:yge
        for k=1:zge
            
            i
            j
            k
            
            
            %curl
            
            %Jr
            
            % xgj(i)
            
            r=xgj(i);
            theta=ygj(j);
            phi=zgj(k);
            

            
            
            %%
            % A diff
            
            Apt = diff_A3d(Bp,2,ygj,[i,j,k]);
            Atp = diff_A3d(Bt,3,zgj,[i,j,k]);
            Arp = diff_A3d(Br,3,zgj,[i,j,k]);
            Apr = diff_A3d(Bp,1,xgj,[i,j,k]);
            Atr = diff_A3d(Bt,1,xgj,[i,j,k]);
            Art = diff_A3d(Br,2,ygj,[i,j,k]);
            
            % Apt= (Bp(i,j+1,k)- Bp(i,j-1,k))/(ygj(j+1)-ygj(j-1));
            % Atp= (Bt(i,j,k+1)- Bt(i,j,k-1))/(zgj(k+1)-zgj(k-1));
            % Arp= (Br(i,j,k+1)- Br(i,j,k-1))/(zgj(k+1)-zgj(k-1));
            % Apr= (Bp(i+1,j,k)- Bp(i-1,j,k))/(xgj(i+1)-xgj(i-1));
            % Atr= (Bt(i+1,j,k)- Bt(i-1,j,k))/(xgj(i+1)-xgj(i-1));
            % Art= (Br(i,j+1,k)- Br(i,j-1,k))/(ygj(j+1)-ygj(j-1));
            
            % seperated differences
            
            Jr(i,j,k)= 1/(r*sin(theta))*(Apt*sin(theta)+Bp(i,j,k)* cos(theta)) - 1/(r*sin(theta))*Atp;
            
            Jt(i,j,k)= 1/(r*sin(theta)) * Arp - 1/r* (Bp(i,j,k)+r*Apr)  ;
            
            Jp(i,j,k)= 1/r* (Bt(i,j,k)+r*Atr) - 1/r* Art;
            
            
            [Jxyz]=datas2c(Jr(i,j,k),Jt(i,j,k),Jp(i,j,k),pi/2-theta, phi);
            
            Jx(i,j,k)=Jxyz(1);
            Jy(i,j,k)=Jxyz(2);
            Jz(i,j,k)=Jxyz(3);
            
            %% write hear
            
        end
        
    end
end

mu0=4*pi*10^-7;

curlr=Jr/mu0*10^9;
curlt=Jt/mu0*10^9;
curlp=Jp/mu0*10^9;

curlx=Jx/mu0*10^9;
curly=Jy/mu0*10^9;
curlz=Jz/mu0*10^9;


e=1.6*10^-19;
me=9.1*10^-31;
mi=1.67*10^-27; 
Nco2=[1*10^10 ] *10^6;
Ne=[1*10^5 ] *10^6; %n_O2+


%load Ndata.mat
%load Ndata2023.mat

%Ne1=squeeze(QuanN3(1,:,:));
%Ne2=[Ne1 ] *10^6;
%Ne2=( [QuanN3(1,:,:) ] )*10^6;

%Me=mean(mean(Ne2));

%Me1=( 2500 )*10^6;


%windr=Jr/mu0./Me1/e;
%windt=Jt/mu0./Me1/e;
%windp=Jp/mu0./Me1/e;

%da=1;
%k = (1/da^2)*ones(da);


%% Tu iono
%%

figure
pos=[1 1 30 16];

set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,1,[.1 .1],[.05 .05],[.06 .1])

 axes(ha(1))
    
    ax1=gca;


qx=curlr(1,:,:);
qr=curlr(1,:,:);
qt=curlt(1,:,:);
qp=curlp(1,:,:);

%qx=windr(1,:,:);
%qr=windr(1,:,:);
%qt=windt(1,:,:);
%qp=windp(1,:,:);



ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on','MapLatLimit',[-85 85])
  
  axesm('MapProjection','eqdcylin','Frame','on','MapLatLimit',[-85 85])
  
 % axesm('MapProjection','mollweid','Frame','on','MapLatLimit',[-85 85])
  
 % 'MapLonLimit',[-90 90]
  
setm(gca,'FLineWidth',4)

    axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;

p1=plabel;
%{
for i=1:length(p1)
 
   if i<6
  p1(i).Position=p1(i).Position-[0,0.1,0];
   end
   
   if i>6
   p1(i).Position=p1(i).Position+[0,0.1,0];
   end
    p1(i).FontSize=15;
end
%}

p1=plabel('FontSize',15);

m1=mlabel;

for i=1:length(m1)
    m1(i).String{1}=2*i-2;
    m1(i).FontSize=15;
    m1(i).Position=m1(i).Position-[0,3.3,0];
end
    

xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);


mt1=nanmean(qt);

mt1(1:5)=0;
mt1(15:20)=0;

on1=ones(1,10);
pj1=on1'*mt1;

%qt=qt-pj1; % decrease Jcf


%qt=qt-10;

qtot=(qr.^2+qt.^2+qp.^2).^0.5;

% add addtional line
xx1j=[xx1,zeros(sz1(2),1)];
xx1j=[xx1j;zeros(1,sz1(3)+1)];

xx2j=[qtot,zeros(sz1(2),1)];
xx2j=[xx2j;zeros(1,sz1(3)+1)];

[gr1,gr2]=meshgrid(ygj/pi*180,zgj/pi*180);

%q1=quiverm(-70,-60,10,0,'r'); hold on;
%q1=quiverm(-80,-60,10,-30,'g'); hold on;


%pcolorm(Ygrid/pi*180,Zgrid/pi*180,xx1j) ;  hold on;

%hq=quiver3(Xq,Yq,Zq,qx1*0,qy1,qz1,6,'k','LineWidth',1,'MaxHeadSize',0.01); hold on;

da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

nm=(qt1.^2+qp1.^2).^0.5;


sst=conv2(xx2j,k,'same');

xx3j=[-pj1,zeros(sz1(2),1)];
xx3j=[xx3j;zeros(1,sz1(3)+1)];


[c,h]=contourm(Ygrid/pi*180,Zgrid/pi*180,double(sst),0:50:400, 'linewidth',2) ;hold on;


c1=qt1'./nm';
c2=qp1'./nm';


q1=quiverm(gr1,gr2,-qt1',qp1',1); hold on;

%q1=quiverm(gr1,gr2,-pj1',qp1'.*0,1); hold on;

%quiverC2D(gr2,gr1,qp1',-qt1',1); hold on;


%[c,h]=contourm(Ygrid/pi*180,Zgrid/pi*180,sst,'LevelStep',15,'ShowText','on') ;hold on;
%[c,h]=contourm(Ygrid/pi*180,Zgrid/pi*180,double(sst),'ShowText','on','LineWidth',2) ;hold on;

%t = clabelm(c,h);
%set(t,'BackgroundColor','white')
%set(t,'BackgroundColor','none')

q1(1).Color=[0 0 0];
q1(2).Color=[0 0 0];
%q1(1).LineWidth=1;
q1(1).LineWidth=1;
q1(2).LineWidth=1;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx1j,-20:5:20)


ter1=90*ones(1,100);
ter2=-90*ones(1,100);
ter3=0*ones(1,100);

tlat=linspace(-90,90,100);

plotm(tlat,ter1,'Linestyle','--','Color','k','Linewidth',1.5)
plotm(tlat,ter2,'Linestyle','--','Color','k','Linewidth',1.5)
plotm(tlat,ter3,'Linestyle','--','Color','k','Linewidth',1.5)

%title('Wind Speed','FontSize',20)

text(-0.4,-1.75,'Local time (h)','FontSize',20);
text(-3.3,-0.4,'Latitude (\circ)','FontSize',20,'Rotation',90);


%text(-3.3,1.7,'a','FontSize',25,'FontWeight','bold');
%text(-3.3,-2.1,'c','FontSize',25,'FontWeight','bold');


caxis=[0 400];

h1=colorbar('eastoutside');
h1.Position=[0.88 0.2 0.02 0.60]
h1.Label.String = '\it{J}\rm{  (nA m^{-2})}';
%h1.Label.String = '(m/s)';

h1.FontSize=15;
h1.FontWeight='normal';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;


h1.Visible='off'

%colormap(mycolormap9)

%set(gcf,'color','w');
colormap jet



caxis=[0 400];

%h1=colorbar('eastoutside');
%h1.Position=[0.88 0.2 0.02 0.60]
%h1.Label.String = '\it{J}\rm{  (nA m^{-2})}';
%h1.Label.String = '(m/s)';

%h1.FontSize=15;
%h1.FontWeight='normal';
%h1.TickDirection='out';
%h1.TickLength=0.01;
%h1.LineWidth=1;

%colormap(mycolormap9)

%[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,'LevelStep',20,'ShowText','off','LineWidth',2) ; hold on;

hl=clegendm(c,h,1);



%c1.LineWidth=2;
set(gca,'FontSize',16)
%set(gcf,'color','w');
colormap jet

hl.Position=[0.9 0.5512 0.0770 0.3249]
%%

print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig5/Fig5_MSO.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig5/Fig5_MSO.eps']);

