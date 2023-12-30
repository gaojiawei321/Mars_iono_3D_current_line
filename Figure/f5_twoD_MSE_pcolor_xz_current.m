clear all;close all;clc


% pcolor is not right in one grid

%%
%readme=' quanx Bx   quanx By  quanz Bz  '
%save threeDMSO_B_xyz_Bt10 quanx quany quanz readme


%load threeDMSE_B_xyz_Bt10_20221212
load threeDMSE_B_xyz_Bt10_20230718


%load threeDMSE_B_xyz_Bt10_20230711

load mycolormap9
load myredblue

Bx=permute(quanx,[2 1 3]);
By=permute(quany,[2 1 3]);
Bz=permute(quanz,[2 1 3]);

%xgrid=linspace(-2,2,21);

%xjg=[-1.9:0.2:1.9];
xj=Xgrid(2)-Xgrid(1);


xjg=Xgrid(1:end-1)+xj/2;
yjg=Ygrid(1:end-1)+xj/2;
zjg=Zgrid(1:end-1)+xj/2;

[XX,YY,ZZ]=meshgrid(Xgrid,Xgrid,Xgrid);

[X,Y,Z]=meshgrid(xjg,xjg,xjg);
%change RE to M

X1=X.*3393.5*1000;
Y1=Y.*3393.5*1000;
Z1=Z.*3393.5*1000;

mu0=4*pi*10^-7;

[curlx,curly,curlz,cav] = curl(X1,Y1,Z1,Bx.*10^-9,By.*10^-9,Bz.*10^-9)  ;

%change J to J/mu0
% change J (A) to J (nA)
curlx=curlx/mu0*10^9;
curly=curly/mu0*10^9;
curlz=curlz/mu0*10^9;


Curlx=zeros(length(Xgrid),length(Xgrid),length(Xgrid));
Curly=zeros(length(Xgrid),length(Xgrid),length(Xgrid));
Curlz=zeros(length(Xgrid),length(Xgrid),length(Xgrid));

Curlx(1:length(Xgrid)-1,1:length(Xgrid)-1,1:length(Xgrid)-1)=curlx;
Curly(1:length(Xgrid)-1,1:length(Xgrid)-1,1:length(Xgrid)-1)=curly;
Curlz(1:length(Xgrid)-1,1:length(Xgrid)-1,1:length(Xgrid)-1)=curlz;


%J (mA/km^2)
%J (nA/m^2)



%%

%j1=curlx(10,:,:);

%zx1=find(isnan(quanx));
%std(zx1)=0;


%%
figure
pos=[1 1 42 12];

v1=-9;
v2=12;

set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,3,[.15 .05],[.12 .08],[.05 .03]);

%a=[0 1 1; 0 1 0];
%pcolor([1,2,3],[1,2],a);
%% slice 
    axes(ha(1))
    ax1=gca;
    
    dd1=16;
    %dd1=32;
    

    LString= 'J_x (nA/m^2)';
%xjg1=[Xgrid(1)-xj Xgrid];




plot_2d_curr_global_quiv(Curlx,dd1,Xgrid,Ygrid,Zgrid,LString,Curlx,Curlz)
title('J_x')

plot_2d_curr_iono_quiv_xz(1)



%%
    axes(ha(2))
    ax1=gca;
    
  %  dd1=14+1;
    LString= 'J_y (nA/m^2)';
%xjg1=[Xgrid(1)-xj Xgrid];
plot_2d_curr_global_quiv(Curly,dd1,Xgrid,Ygrid,Zgrid,LString,Curlx,Curlz)
title('J_y')

plot_2d_curr_iono_quiv_xz(2)

%%
    axes(ha(3))
    ax1=gca;
    
  %  dd1=14+1;
    LString= 'J_z (nA/m^2)';
%xjg1=[Xgrid(1)-xj Xgrid];
plot_2d_curr_global_quiv(Curlz,dd1,Xgrid,Ygrid,Zgrid,LString,Curlx,Curlz)
title('J_z')

plot_2d_curr_iono_quiv_xz(3)
%colormap(mycolormap9)
colormap(myredblue)


%%
%%
%{
    axes(ha(1))
    ax1=gca;

    xjg=[-3.1:0.2:2.9];
 [X,Y,Z]=meshgrid(xjg,xjg,xjg);

s1=slice(X,Y,Z, Curlx, [],[0],[]); hold on;

%middle of 11 and 12
dd1=15+1;
%px=(Curlx(10,:,:)+Curlx(11,:,:))/2;
%py=(Curly(10,:,:)+Curly(11,:,:))/2;
%pz=(Curlz(10,:,:)+Curlz(11,:,:))/2;

%q1=quiver3(X(dd1,:,:),Y(dd1,:,:),Z(dd1,:,:),Curlx(dd1,:,:),Curly(dd1,:,:),Curlz(dd1,:,:),2,'k'); hold on;

q1=quiver3((X(dd1,:,:)+X(dd1+1,:,:))/2,(Y(dd1,:,:)+Y(dd1+1,:,:))/2,(Z(dd1,:,:)+Z(dd1+1,:,:))/2,(Curlx(dd1,:,:)+Curlx(dd1+1,:,:))/2,(Curly(dd1,:,:)+Curly(dd1+1,:,:))/2*0,(Curlz(dd1,:,:)+Curlz(dd1+1,:,:))/2,2,'k'); hold on;

%q1=quiver((X(dd1,:,:)+X(dd1+1,:,:))/2,(Y(dd1,:,:)+Y(dd1+1,:,:))/2,(Curlx(dd1,:,:)+Curlx(dd1+1,:,:))/2,(Curly(dd1,:,:)+Curly(dd1+1,:,:))/2,2,'k'); hold on;

q1=quiver3((X(dd1,:,:)+X(dd1+1,:,:))/2,(Y(dd1,:,:)+Y(dd1+1,:,:))/2,(Z(dd1,:,:)+Z(dd1+1,:,:))/2,(Curlx(dd1,:,:)+Curlx(dd1+1,:,:))/2,(Curly(dd1,:,:)+Curly(dd1+1,:,:))/2,(Curlz(dd1,:,:)+Curlz(dd1+1,:,:))/2,2,'k'); hold on;


%[startX,startY,startZ] = meshgrid(-1.5,-2:0.2:2,-2:0.2:2);
%[startX,startY,startZ] = meshgrid([-1,0,1],0,-2:0.2:2);

%[startX,startZ] = meshgrid([-1,0,1],-2:0.2:2);

%[X1,Z1]=meshgrid(xgrid,xgrid);
%verts = stream2(X1,Z1,squeeze(curlx(:,10,:)),squeeze(curlz(:,10,:)),startX,startZ);

%verts = stream3(X,Y,Z,curlx,curly,curlz,startX,startY,startZ);

%s1=streamline(verts); hold on;

%for ci=1:length(s1)
%    s1(ci).Color=[0 0 0];
%end

%quiver3(X1,Y1,Z1,Bx,By,Bz); hold on;

%quiver3(X,Y,Z,curlx,curly,curlz,'k');

for kj=1:length(xslice)
    xkj=xslice(kj);
plot3([xkj,xkj,xkj,xkj,xkj],[-2,2,2,-2,-2],[-2,-2,2,2,-2],'k','LineWidth',1) ;hold on;
end


xlabel('X_{MSE} (R_m)')
ylabel('Y_{MSE} (R_m)')
zlabel('Z_{MSE} (R_m)')
caxis([-10,10]);
h1=colorbar;
h1.Label.String = 'J_x (nA/m^2)';
 shading flat
 title('Bcrust<10 nT  J_x')
% title('Bcrust<10 nT  B_{\theta} ')
% title('Bcrust<10 nT  B_{\rho} std')

 
 set(gca,'YDir','reverse');
 set(gca,'XDir','reverse');

 %view(76,13)
 % view(v1,v2)
  
  view(0,0)
  
xlim([-1.5 1.5])
ylim([-1.8 1.8])
zlim([-1.8 1.8])

 %colormap(mycolormap7)
 set(gca,'linewidth',1.5)
 
set(gca,'Fontsize',10);
 

%%

%%

    axes(ha(2))
    ax1=gca;

%slice(X,Y,Z, curly, xslice,[],[]); hold on;
slice(X,Y,Z, Curly, [],[0],[]); hold on;

q1=quiver3(X(dd1,:,:),Y(dd1,:,:),Z(dd1,:,:),Curlx(dd1,:,:),Curly(dd1,:,:),Curlz(dd1,:,:),2,'k'); hold on;

%q1=quiver3(X(dd1,:,:)+0.1,Y(dd1,:,:)+0.1,Z(dd1,:,:)+0.1,curlx(dd1,:,:),curly(dd1,:,:),curlz(dd1,:,:),2,'k'); hold on;

for kj=1:length(xslice)
    xkj=xslice(kj);
plot3([xkj,xkj,xkj,xkj,xkj],[-2,2,2,-2,-2],[-2,-2,2,2,-2],'k','LineWidth',1) ;hold on;
end


xlabel('X_{MSE} (R_m)')
ylabel('Y_{MSE} (R_m)')
zlabel('Z_{MSE} (R_m)')

caxis([-10,10]);
h1=colorbar;
h1.Label.String = 'J_y (nA/m^2)';
 shading flat
 title('Bcrust<10 nT  J_y')
% title('Bcrust<10 nT  B_{\theta} ')
% title('Bcrust<10 nT  B_{\rho} std')

 
 set(gca,'YDir','reverse');
 set(gca,'XDir','reverse');

 %view(76,13)
  view(v1,v2)
  
    view(0,0)
  
xlim([-1.5 1.5])
ylim([-1.8 1.8])
zlim([-1.8 1.8])

 %colormap(mycolormap7)
 set(gca,'linewidth',1.5)
set(gca,'Fontsize',10);
 a=1
 
 %%
 
    axes(ha(3))
    ax1=gca;

%slice(X,Y,Z,curlz, xslice,[],[]); hold on;
slice(X,Y,Z, Curlz, [],[0],[]); hold on;

%q1=quiver3(X(dd1,:,:)+0.1,Y(dd1,:,:)+0.1,Z(dd1,:,:)+0.1,Curlx(dd1,:,:),Curly(dd1,:,:),Curlz(dd1,:,:),2,'k'); hold on;
q1=quiver3(X(dd1,:,:),Y(dd1,:,:),Z(dd1,:,:),Curlx(dd1,:,:),Curly(dd1,:,:),Curlz(dd1,:,:),2,'k'); hold on;

%quiver3(X,Y,Z,Bx,By,Bz);


for kj=1:length(xslice)
    xkj=xslice(kj);
plot3([xkj,xkj,xkj,xkj,xkj],[-2,2,2,-2,-2],[-2,-2,2,2,-2],'k','LineWidth',1) ;hold on;
end


xlabel('X_{MSE} (R_m)')
ylabel('Y_{MSE} (R_m)')
zlabel('Z_{MSE} (R_m)')

caxis([-10,10]);
h1=colorbar;
h1.Label.String = 'J_z (nA/m^2)';
 shading flat
 title('Bcrust<10 nT  J_z')
% title('Bcrust<10 nT  B_{\theta} ')
% title('Bcrust<10 nT  B_{\rho} std')

 
 set(gca,'YDir','reverse');
 set(gca,'XDir','reverse');

 %view(76,13)
 % view(v1,v2)
    view(0,0)
  
xlim([-1.5 1.5])
ylim([-1.8 1.8])
zlim([-1.8 1.8])

 %colormap(mycolormap7)
 set(gca,'linewidth',1.5)
set(gca,'Fontsize',10);
 
 a=1
  
 colormap(mycolormap8)
 
 %}
 %%
 
%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Bx_slice_3.tiff']);
%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Btheta_slice_3.tiff']);
%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Brho_slice_3_std.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig4/Jmse_2d_y0_3.eps']);
print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig4/Jmse_2d_y0_3.tiff']);


print('-depsc','-r300',['C:\Users\Gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Jmse_2d_y0_3.eps']);
print('-dtiff','-r300',['C:\Users\Gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Jmse_2d_y0_3.tiff']);



 