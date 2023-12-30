clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate current in Martian ionosphere
% 2022_10_30 modified JW Gao


%load ionMSEsph_Brtp_L10_-F50_20221113
%load ionMSEsph_Brtp_Bt10_20221109

%load ionMSEsph_Brtp_L10_20221113
%load ionMSEsph_Brtp_L10_Bm150_20221220

%load ionMSEsph_Brtp_L20_Bm120_20230312
%load ionMSEsph_Brtp_L10_Bm120_Amse_20230312
load ionMSOsph_Brtp_L10_Bm120_2023

%load ionMSEsph_Brtp_L10_Bm120_Amse_20221222


%load ionMSEsph_Brtp_L10_Bm150_20221220

load mycolormap9
%%
method='linear';
quanBx=fillmissing(quanBx,method);
quanBy=fillmissing(quanBy,method);
quanBz=fillmissing(quanBz,method);


quanBr=fillmissing(quanBr,method);
quanBt=fillmissing(quanBt,method);
quanBp=fillmissing(quanBp,method);

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


%% using curl not good
%{
%[X,Y,Z]=meshgrid(xgj,ygj,zgj);

%[x,y,z]=sph2cart(Z,Y,X);

%[curlx,curly,curlz,cav] = curl(x,y,z,quanBx.*10^-9,quanBy.*10^-9,quanBz.*10^-9)  ;

for i=1:xge
    for j=1:yge
        for k=1:zge

        [x,y,z]=sph2cart(zgj(k),ygj(j),xgj(i));

        Cartx(i,j,k)=x;
        Carty(i,j,k)=y;
        Cartz(i,j,k)=z;

        end
    end
end

[curlx,curly,curlz,cav] = curl(Cartx,Carty,Cartz,quanBx.*10^-9,quanBy.*10^-9,quanBz.*10^-9)  ;
%}


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
            
            % seems same? double-sided differences
            %  d1=(Bp(i,j+1,k)- Bp(i,j,k))/(ygj(j+1)-ygj(j))*sin(theta) + Bp(i,j,k)* cos(theta);
            
            
            %   d1=(Bp(i,j+1,k)*sin(ygj(j+1)) - Bp(i,j-1,k)*sin(ygj(j-1)) )/(ygj(j+1)-ygj(j-1)) ;
            %   d2=(Bt(i,j,k+1)-Bt(i,j,k-1))/(zgj(k+1)-zgj(k-1));
            %   d3=(Br(i,j,k+1)-Br(i,j,k-1))/(zgj(k+1)-zgj(k-1));
            %   d4=(Bp(i+1,j,k)*xgj(i+1) - Bp(i-1,j,k)*xgj(i-1) )/(xgj(i+1)-xgj(i-1)) ;
            %   d5=(Bt(i+1,j,k)*xgj(i+1) - Bt(i-1,j,k)*xgj(i-1) )/(xgj(i+1)-xgj(i-1)) ;
            %   d6=(Br(i,j+1,k)-Br(i,j-1,k))/(ygj(i+1)-ygj(i-1)) ;
            
            
            %  diff single-sided differences
            %  d1=(Bp(i,j+1,k)*sin(ygj(j+1)) - Bp(i,j,k)*sin(ygj(j)))/(ygj(j+1)-ygj(j)) ;
            %  d2=(Bt(i,j,k+1)-Bt(i,j,k))/(zgj(k+1)-zgj(k));
            %  d3=(Br(i,j,k+1)-Br(i,j,k))/(zgj(k+1)-zgj(k));
            %  d4=(Bp(i+1,j,k)*xgj(i+1) - Bp(i,j,k)*xgj(i) )/(xgj(i+1)-xgj(i)) ;
            %  d5=(Bt(i+1,j,k)*xgj(i+1) - Bt(i,j,k)*xgj(i) )/(xgj(i+1)-xgj(i)) ;
            %  d6=(Br(i,j+1,k)-Br(i,j,k))/(ygj(i+1)-ygj(i)) ;
            
            
            %  Jr(i,j,k)= 1/(r*sin(theta))*d1 - 1/(r*sin(theta))*d2;
            
            %  Jt(i,j,k)= 1/(r*sin(theta)) * d3 - 1/r* d4  ;
            
            %  Jp(i,j,k)= 1/r* d5 - 1/r* d6;
            
            
            
            
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
Ne=[1*10^4 ] *10^6;

%windr=Jr/mu0/Ne/e;
%windt=Jt/mu0/Ne/e;
%windp=Jp/mu0/Ne/e;

windr=curlr;
windt=curlt;
windp=curlp;



%% Tu iono figure 5b
jiang=0:50:400;

figure
pos=[1 1 40 20];

set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,2,[.05 .05],[.05 .05],[.05 .05])

 axes(ha(1))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])

%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[-90 90]) 

%axesm('MapProjection','ortho','Origin',[0,0,0],'Grid','on','Frame','on','MlineLimit',[-75 75]);

%axesm('MapProjection','ortho','Origin',[0,90,0],'Grid','on','Frame','on','MlineLimit',[60 90]);

axesm('ortho','maplatlim',[50 90]);

%axesm('eqaazim','MapLatLimit',[0 90])

  
  
setm(gca,'FLineWidth',4)

 %   axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;
plabel on;
mlabel on;
setm(gca,'MLabelParallel',0,'PLabelMeridian',0)

xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

q1=quiverm(gr1,gr2,-qt1',qp1',2,'AutoScaleFactor',0.9,'MaxHeadSize',0.2); hold on;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

%[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,'LevelStep',20,'ShowText','off','LineWidth',2) ; hold on;

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;

%'LevelStep',10

%t = clabelm(c,h);
%c1=clegendm(c,h,0);
%c1.LineWidth=0.5;
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
set(gca,'FontSize',16)

%LevelStep',40

%%

 axes(ha(2))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])
%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[90 -90]) 
axesm('ortho','maplatlim',[-90 -50]);

%axesm('MapProjection','ortho','Origin',[0,0,0],'Grid','on','Frame','on','MlineLimit',[-75 75]);


%axesm('eqaazim','MapLatLimit',[0 90])

  
  
setm(gca,'FLineWidth',4)

 %   axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;
plabel on;
mlabel on;
setm(gca,'MLabelParallel',0,'PLabelMeridian',0)

xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

%q1=quiverm(gr1,gr2,-qt1',qp1',2); hold on;
q1=quiverm(gr1,gr2,-qt1',qp1',2,'AutoScaleFactor',0.9,'MaxHeadSize',0.2); hold on;


%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

%[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,'LevelStep',20,'ShowText','off','LineWidth',2) ; hold on;

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;
%[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,'LevelStep',5,'ShowText','off','LineWidth',1.5) ; hold on;

%t = clabelm(c,h);
%clegendm(c,h,0)
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
%LevelStep',40

set(gca,'FontSize',16)

%set(gca,'Linewidth',2)

colormap(jet)

%print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5bz.tiff']);
%print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5bz.eps']);

print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig5/Fig5b_MSO.tiff']);
print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig5/Fig5b_MSO.eps']);



%% figure 3D sphere



figure
pos=[1 1 40 20];

set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,2,[.05 .05],[.05 .05],[.05 .05])

 axes(ha(1))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])
%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[-90 90]) 
axesm('MapProjection','ortho','Origin',[45,45,-30],'Grid','on','Frame','on','MlineLimit',[-75 75]);


%axesm('eqaazim','MapLatLimit',[0 90])

setm(gca,'FLineWidth',4)

    axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;


h1=plabel('FontSize',15,'FontWeight','bold');
h2=mlabel('FontSize',15,'FontWeight','bold');

%m1=mlabel;
%p1=plabel;

%for i=1:length(m1)
%    m1(i).FontWeight='bold';
%    m1(i).FontSize=18;
%end
    
%for i=1:length(p1)
%    p1(i).FontWeight='bold';
%    p1(i).FontSize=18;
%end
    


setm(gca,'MLabelParallel',0,'PLabelMeridian',0)


xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

q1=quiverm(gr1,gr2,-qt1',qp1',1.5); hold on;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;

%'LevelStep',10

%t = clabelm(c,h);
c1=clegendm(c,h,-1);
%c1.LineWidth=0.5;
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
set(gca,'FontSize',16)

%LevelStep',40


 axes(ha(2))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])
%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[90 -90]) 

axesm('MapProjection','ortho','Origin',[-45,-135,30],'Grid','on','Frame','on','MlineLimit',[-75 75]);

%axesm('eqaazim','MapLatLimit',[0 90])
    axis off
  
  
setm(gca,'FLineWidth',4)

 %   axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;
h1=plabel('FontSize',15,'FontWeight','bold');
h2=mlabel('FontSize',15,'FontWeight','bold');

setm(gca,'MLabelParallel',0,'PLabelMeridian',0)

xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

q1=quiverm(gr1,gr2,-qt1',qp1',1.5); hold on;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;
%t = clabelm(c,h);
clegendm(c,h,-1)
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
%LevelStep',40

set(gca,'FontSize',16)

%set(gca,'Linewidth',2)

colormap(jet)
%%

print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5b3.tiff']);

print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5b3.eps']);


%%
%% figure 3D sphere


figure
pos=[1 1 40 20];

set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,2,[.05 .05],[.05 .05],[.05 .05])

 axes(ha(1))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])
%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[-90 90]) 
axesm('MapProjection','ortho','Origin',[-45,45,30],'Grid','on','Frame','on','MlineLimit',[-75 75]);


%axesm('eqaazim','MapLatLimit',[0 90])

setm(gca,'FLineWidth',4)

    axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;


h1=plabel('FontSize',15,'FontWeight','bold');
h2=mlabel('FontSize',15,'FontWeight','bold');

%m1=mlabel;
%p1=plabel;

%for i=1:length(m1)
%    m1(i).FontWeight='bold';
%    m1(i).FontSize=18;
%end
    
%for i=1:length(p1)
%    p1(i).FontWeight='bold';
%    p1(i).FontSize=18;
%end
    


setm(gca,'MLabelParallel',0,'PLabelMeridian',0)


xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

q1=quiverm(gr1,gr2,-qt1',qp1',1.5); hold on;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;

%'LevelStep',10

%t = clabelm(c,h);
c1=clegendm(c,h,-1);
%c1.LineWidth=0.5;
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
set(gca,'FontSize',16)

%LevelStep',40


 axes(ha(2))
    
    ax1=gca;


qx=windr(1,:,:);

qr=windr(1,:,:);
qt=windt(1,:,:);
qp=windp(1,:,:);


ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

sz1=size(qx);
  %  axesm mercator ;

%  axesm('MapProjection','miller','Frame','on')

%axesm('eqdazim','MapLonLimit',[-90 90])
%axesm('eqdazim','FLatLimit',[],'MapLonLimit',[90 -90]) 

axesm('MapProjection','ortho','Origin',[45,-135,-30],'Grid','on','Frame','on','MlineLimit',[-75 75]);

%axesm('eqaazim','MapLatLimit',[0 90])
    axis off
  
  
setm(gca,'FLineWidth',4)

 %   axis off
    dc=0;
framem on; %gridm
gridm on;
%pi=plabel;
h1=plabel('FontSize',15,'FontWeight','bold');
h2=mlabel('FontSize',15,'FontWeight','bold');

setm(gca,'MLabelParallel',0,'PLabelMeridian',0)

xx1 =squeeze(qr);
qr=squeeze(qr);
qt=squeeze(qt);
qp=squeeze(qp);

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


qt(isnan(qt))=0;
qp(isnan(qp))=0;


da=2;
k = (1/da^2)*ones(da);
qt1=conv2(qt,k,'same');
qp1=conv2(qp,k,'same');

q1=quiverm(gr1,gr2,-qt1',qp1',1.5); hold on;

%contourm(Ygrid/pi*180,Zgrid/pi*180,xx2j) ;hold on;

q1(1).Color=[0 0 0 ];
q1(2).Color=[0 0 0 ];
%q1(1).LineWidth=1;
q1(1).LineWidth=2;
q1(2).LineWidth=2;

%k=1/9*ones(3);
k = (1/4)*ones(2);

sst=conv2(xx2j,k,'same');

[c,h] = contourm(Ygrid/pi*180,Zgrid/pi*180,sst,jiang,'ShowText','off','LineWidth',1.5) ; hold on;
%t = clabelm(c,h);
clegendm(c,h,-1)
%set(t,'FontWeight','bold')
%set(t,'BackgroundColor','none')
%LevelStep',40

set(gca,'FontSize',16)

%set(gca,'Linewidth',2)

colormap(jet)
%%

print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5b4.tiff']);

print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig5\Fig5b4.eps']);


%%