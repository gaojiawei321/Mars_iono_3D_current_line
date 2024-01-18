clear all;close all;clc

%%%%%%%%%% model at 150 km

rm=3393.5;
aj=1;
ej=1;
Az=aj:aj:360;
El=-90+ej:ej:90;
az=Az/180*pi;
el=El/180*pi;

[AZ,EL]=meshgrid(az,el);

al=120;
h=rm+al*ones(length(el),length(az));

[lx,ly,lz] = sph2cart(AZ,EL,h);

Lx=reshape(lx,1,[]);
Ly=reshape(ly,1,[]);
Lz=reshape(lz,1,[]);

a=[Lx;Ly;Lz];

brtp=mgao_r(a,'g_110_mm_q','h_110_mm_q');
%brtp=gao_r_c(a);

%brtp=lang18_r(a);
%brtp=cain_r(a);



for i=1:3
  %  bx=reshape(bca(1,:),length(el),[]);
  %  by=reshape(bca(2,:),length(el),[]);
  %  bz=reshape(bca(3,:),length(el),[]);

    br=reshape(brtp(1,:),length(el),[]);
    bt=reshape(brtp(2,:),length(el),[]);
    bp=reshape(brtp(3,:),length(el),[]);
    
    
end

ball=(br.^2+bt.^2+bp.^2).^0.5;

bl1=find(ball<10);

%ball(bl1);



a=1


%%




%%    



%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate current in Martian ionosphere
% 2022_11_13

%
%load MAVEN_4s_q_2022.mat

load MAVEN_4s_q_202211.mat

%load MAVEN_4s_f.mat

Rm=3393.5;


%[theta,rho,z] = cart2pol(MAVEN_4s_f.Lmse_zh(:,3),MAVEN_4s_f.Lmse_zh(:,2),MAVEN_4s_f.Lmse_zh(:,1));
%bpol=datac2p(MAVEN_4s_f.Bmse_zh(:,3),MAVEN_4s_f.Bmse_zh(:,2),MAVEN_4s_f.Bmse_zh(:,1),theta')';


Bmodel=MAVEN_4s_f.Bpc-MAVEN_4s_f.Bpc_res_g;
%Bmodel=MAVEN_4s_f.Bpc-MAVEN_4s_f.Bpc_res_l18;

Bot=vecnorm(MAVEN_4s_f.Bmso',2);


Bmt=(Bmodel(:,1).^2+Bmodel(:,2).^2+Bmodel(:,3).^2).^0.5;
%Bimf=vecnorm(MAVEN_4s_f.IMFB')';
%Vimf=vecnorm(MAVEN_4s_f.IMFV')';


[pcp, pct, pcr]=cart2sph(MAVEN_4s_f.Lpc(:,1),MAVEN_4s_f.Lpc(:,2),MAVEN_4s_f.Lpc(:,3));

zz1=find(pcp<0);
azimuth=pcp;
azimuth(zz1)=2*pi+pcp(zz1);

azimuth=azimuth/pi*180;
azimuth=round(azimuth);
zz2=find(azimuth==360);
azimuth(zz2)=0;
elevation=round(pct/pi*180);




%Bmodel150=ball(elevation(1:100)+90,azimuth(1:100)+1);
Bmodel150=zeros(length(elevation),1);
for i=1:length(elevation)
    Bmodel150(i)=ball(elevation(i)+90,azimuth(i)+1);
end

    

L=( MAVEN_4s_f.Lmso(:,1).^2+MAVEN_4s_f.Lmso(:,2).^2+MAVEN_4s_f.Lmso(:,3).^2).^0.5-Rm;

l1=0;
l2=600;
%hou=find(Bmt<10);


%hou=find(Bmt<10& L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 );

%hou=find(Bmt<10& L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 &MAVEN_4s_f.mse_pan==1 );

%hou=find(Bmt<10& L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 &MAVEN_4s_f.mse_pan==1 & 5.*Bmt<Bot');
%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 &MAVEN_4s_f.mse_pan==1 & 5.*Bmt<Bot');


%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 &MAVEN_4s_f.mse_pan==1);

%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 & MAVEN_4s_f.mse_pan==1);

hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 & MAVEN_4s_f.mse_pan==1 & MAVEN_4s_f.F107>50);

% add z hemisphere 
%hou=find(Bmt<10 & L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 & MAVEN_4s_f.Lpc(:,3)>0 );

[La,Le,Lr] = cart2sph( MAVEN_4s_f.Lmse_zh(hou,1),MAVEN_4s_f.Lmse_zh(hou,2),MAVEN_4s_f.Lmse_zh(hou,3));

Let=pi/2-Le;

%B in MSE; cart to sph 
%% check this on: seem make sence


%MAVEN_4s_f.Bmse_zh(hou(1),:)
%MAVEN_4s_f.Lmse_zh(hou(1),:)
%[La(1),Le(1),Lr(1)]/pi*180
%Bmse_sph1=datac2s(MAVEN_4s_f.Bmse_zh(hou(1),1),MAVEN_4s_f.Bmse_zh(hou(1),2),MAVEN_4s_f.Bmse_zh(hou(1),3),Let(1)',La(1)');
%Bmse_sph1=datac2s(1,0,0,pi/2,pi/2);


%check complete



Bmse_sph=datac2s(MAVEN_4s_f.Bmse_zh(hou,1),MAVEN_4s_f.Bmse_zh(hou,2),MAVEN_4s_f.Bmse_zh(hou,3),Let',La');


a=1

%%%%%% Bmse sph (BR BT BP)   Location mse_sph (LR LT LP)
%%%%%%  Br(ourward) Btheta(southward) Bphi(eastward)

bsph=[ Bmse_sph', Lr, Le, La];

bxyz=[ MAVEN_4s_f.Bmse_zh(hou,:), Lr, Le, La];


clear MAVEN_4s_f
%%
load /Users/gaojiawei/program/mars/crust/colormap_1/mycolormap7


 %%
%z1=find(b(:,1)==0);
%b(z1,:)=[];

xge=10;
yge=10;
zge=20;

xs=100+Rm;
xe=600+Rm;

ys=-pi/2;
ye=pi/2;

zs=-pi;
ze=pi;

xj=(xe-xs)/xge;
yj=(ye-ys)/yge;
zj=(ze-zs)/zge;


%Xgrid=[100  200 300 400 500];

%Xgrid=logspace(log10(xs),log10(xe),xge+1);

%X1=linspace(log10(xs-Rm),log10(xe-Rm),xge+1);

%Xgrid=10.^(X1)+Rm;

%X1=exp(Xgrid)+Rm;

%Xgrid=logspace(xs,xe,xge+1);

Xgrid=linspace(xs,xe,xge+1);
Ygrid=linspace(ys,ye,yge+1);
Zgrid=linspace(zs,ze,zge+1);

xgj=Xgrid(1:end-1)+xj/2;
ygj=Ygrid(1:end-1)+yj/2;
zgj=Zgrid(1:end-1)+zj/2;

quanBr=zeros(xge,yge,zge);
quanBt=zeros(xge,yge,zge);
quanBp=zeros(xge,yge,zge);

quanBx=zeros(xge,yge,zge);
quanBy=zeros(xge,yge,zge);
quanBz=zeros(xge,yge,zge);

ge=zeros(xge,yge,zge);

for i=1:xge
    for j=1:yge
        for k=1:zge
        
        i
        j
        k

               
       z1=find(bsph(:,4)<Xgrid(i+1) &bsph(:,4)>Xgrid(i)    &bsph(:,5)<Ygrid(j+1)&bsph(:,5)> Ygrid(j)  &bsph(:,6)<Zgrid(k+1)  &bsph(:,6)>Zgrid(k) ) ;
         
        %median
        
        quanBr(i,j,k)=nanmean(bsph(z1,1));
        quanBt(i,j,k)=nanmean(bsph(z1,2));
        quanBp(i,j,k)=nanmean(bsph(z1,3));
        
        quanBx(i,j,k)=nanmean(bxyz(z1,1));
        quanBy(i,j,k)=nanmean(bxyz(z1,2));
        quanBz(i,j,k)=nanmean(bxyz(z1,3));
        
        
        
        ge_sph(i,j,k)=length(bsph(z1,3));
        ge_xyz(i,j,k)=length(bxyz(z1,3));
        
        end

    end
end


%%


%Xj=[-1.8:0.2:2];

readme='ion MSE sph quanBr Br   quanBt Bt  quanBp Bp  '

%save /Users/gaojiawei/program/mars/external/quan/tu/ion_mat/ionMSEsph_Brtp_L10_20221113 quanBr quanBt quanBp quanBx quanBy quanBz ge_xyz ge_sph Xgrid Ygrid Zgrid readme

save /Users/gaojiawei/program/mars/external/quan/tu/ion_mat/ionMSEsph_Brtp_L10_+F50_20221113 quanBr quanBt quanBp quanBx quanBy quanBz ge_xyz ge_sph Xgrid Ygrid Zgrid readme

%% 

%}

load ionMSEsph_Brtp_L10_20221113
load /Users/gaojiawei/program/mars/crust/colormap_1/mycolormap8

yge=10;
zge=20;

% PCOLOR
%qz=permute(quanBp,[2 1 3]);

%%
figure
pos=[1 1 35 20];


set(gcf,'unit','centimeters','position',pos)


for hc=1:xge
    
    
ha = tight_subplot(2,2,[.15 .1],[.09 .05],[.08 .05])
%%%

    axes(ha(1))
    ax1=gca;

  %  hc=4;
   % xx1=squeeze(nanmean(quanBx,1));
xx1 =squeeze(quanBx(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_x')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')

    axes(ha(2))
    ax2=gca;

   % xx1=squeeze(nanmean(quanBy,1));
xx1=squeeze(quanBy(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_y')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')




    axes(ha(3))
    ax3=gca;

   % xx1=squeeze(nanmean(quanBz,1));
xx1=squeeze(quanBz(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_z')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')



    axes(ha(4))
    ax4=gca;

    xx1=squeeze(nansum(ge_xyz,1));
   % xx1=squeeze(nanmean(ge,1));
    
%xx=squeeze(quanBx(1,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([0,10000]);
h1=colorbar;
h1.Label.String = '(Counts)';
g1=xgj-Rm;
title(['Altitude=',num2str(g1(hc)),'   (Counts)'])
 

 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')




%%

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/ion/Ion_MSE_xyz_Bt10_1.tiff']);

cun1=['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/ion/B/Ion_MSE_xyz_Bt10_h',num2str(hc),'.tiff']

print('-dtiff ','-r300',cun1);

clf


end



%% second figure


figure
pos=[1 1 35 20];


set(gcf,'unit','centimeters','position',pos)

for hc=1:xge
    
ha = tight_subplot(2,2,[.15 .1],[.09 .05],[.08 .05])
%%%

    axes(ha(1))
    ax1=gca;

 %   xx1=squeeze(nanmean(quanBr,1));
xx1=squeeze(quanBr(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_r')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')

    axes(ha(2))
    ax2=gca;

 %   xx1=squeeze(nanmean(quanBt,1));
xx1=squeeze(quanBt(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_\theta')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')




    axes(ha(3))
    ax3=gca;

%    xx1=squeeze(nanmean(quanBp,1));

xx1=squeeze(quanBp(hc,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([-30,30]);
h1=colorbar;
h1.Label.String = '(nT)';
title('B_\phi')
 

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')



    axes(ha(4))
    ax4=gca;

    xx1=squeeze(nansum(ge_sph,1));
   % xx1=squeeze(nanmean(ge,1));
    
%xx=squeeze(quanBx(1,:,:));

% add addtional line
xx1j=[xx1,zeros(yge,1)];
xx1j=[xx1j;zeros(1,zge+1)];


pcolor(Zgrid,Ygrid,xx1j)
colormap(mycolormap8)
shading flat

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'0','6','12','18','24'})
yticks([-pi/2 -pi/4 0 pi/4 pi/2])
yticklabels({'-90','-45','0','45','90'})

xlabel('Local time (h)')
ylabel('Latitude (°)')

caxis([0,10000]);
h1=colorbar;
h1.Label.String = '(Counts)';


g1=xgj/1000-Rm;
title(['Altitude=',num2str(g1(hc)),'   (Counts)'])

set(gca,'TickDir','out')
    
  set(gca,'linewidth',1.5)
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')




%%

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/ion/Ion_MSE_rtp_Bt10_1.tiff']);

cun1=['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/ion/B/Ion_MSE_rtp_Bt10_h',num2str(hc),'.tiff']


print('-dtiff ','-r300',cun1);

clf

end
 