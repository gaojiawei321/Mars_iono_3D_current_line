clear all;close all;clc


%%
load threeDMSE_B_xyz_Bt10_2023

xgj=Xgrid(1:end-1)+(Xgrid(2)-Xgrid(1))/2;
ygj=Ygrid(1:end-1)+(Ygrid(2)-Ygrid(1))/2;
zgj=Zgrid(1:end-1)+(Zgrid(2)-Zgrid(1))/2;

%load threeDMSE_B_xyz_Bt10
load mycolormap8

qx=permute(quanx,[2 1 3]);
qy=permute(quany,[2 1 3]);
qz=permute(quanz,[2 1 3]);

%%

figure

pos=[1 1 30 20];
%pos=[1 1 30 25];


set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(3,1,[.06 .02],[.08 .05],[.09 .06])


axes(ha(1))
ax1=gca;

plot_3d_mag_global_quiv(qx,Xgrid,Ygrid,Zgrid,'\it{B_x}\rm{  (nT)}',qx,qy,qz)

%plot_3d_mag_global(qx,Xgrid,Ygrid,Zgrid,'B_x (nT)')

text(2.05,-1,4.2,'a','Fontsize',20,'FontWeight','bold','FontName','helvetica')

xlabel([])

axes(ha(2))
ax1=gca;

plot_3d_mag_global_quiv(qy,Xgrid,Ygrid,Zgrid,'\it{B_y}\rm{  (nT)}',qx,qy,qz)

%plot_3d_mag_global(qy,Xgrid,Ygrid,Zgrid,'B_y (nT)')


text(2.05,-1,4.2,'b','Fontsize',20,'FontWeight','bold','FontName','helvetica')

xlabel([])

axes(ha(3))
ax1=gca;

plot_3d_mag_global_quiv(qz,Xgrid,Ygrid,Zgrid,'\it{B_z}\rm{  (nT)}',qx,qy,qz)

xlabel('')
%plot_3d_mag_global(qz,Xgrid,Ygrid,Zgrid,'B_z (nT)')

text(2.05,-1,4.2,'c','Fontsize',20,'FontWeight','bold','FontName','helvetica')

%xlabel([])

load mycolormap12
colormap(mycolormap12)

%%



  
 %%
 
 
 

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Bx_slice_3.tiff']);
%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Btheta_slice_3.tiff']);
%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Brho_slice_3_std.tiff']);

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Bmse_xyz_slice_3_0.5_Vimf<400.tiff']);

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/Bmse_xyz_slice_3_0.5_Bimf>2.tiff']);

%print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig1\Bmse_xyz_slice_3_0.5_Bt10_q1.tiff']);

%print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\Fig1\Bmse_xyz_slice_3_0.5_Bt10_q1.eps']);

print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Bmse_xyz_slice_3_0.5_Bt10_q1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Bmse_xyz_slice_3_0.5_Bt10_q1.eps']);

 
