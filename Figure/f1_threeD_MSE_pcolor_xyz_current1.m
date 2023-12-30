clear all;close all;clc


%%
%readme=' quanx Bx   quanx By  quanz Bz  '
%save threeDMSO_B_xyz_Bt10 quanx quany quanz readme

%load('/Users/gaojiawei/program/mars/external/quan/tu/global_mat/threeDMSE_B_xyz_Bt10_20221107.mat')
load threeDMSE_B_xyz_Bt10_20221212.mat
%load threeDMSE_B_xyz_Bt10_2023.mat



%load('/Users/gaojiawei/program/mars/external/quan/tu/global_mat/threeDMSE_B_xyz_Bt10time_20221112.mat')
%load('/Users/gaojiawei/program/mars/external/quan/tu/global_mat/threeDMSE_B_xyz_Bt10_Con_+120_20221110.mat')
%load('/Users/gaojiawei/program/mars/external/quan/tu/global_mat/threeDMSE_B_xyz_Bt10_Psw_-0.5_20221110.mat')


load mycolormap8

Bx=permute(quanx,[2 1 3]);
By=permute(quany,[2 1 3]);
Bz=permute(quanz,[2 1 3]);

xgrid=[-2.9:0.2:2.9];
[X,Y,Z]=meshgrid(xgrid,xgrid,xgrid);
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

%J (mA/km^2)
%J (nA/m^2)



%%

%j1=curlx(10,:,:);

%zx1=find(isnan(quanx));
%std(zx1)=0;



%save current_MSE_global_1213 curlx curly curlz X1 Y1 Z1

%%



figure
pos=[1 1 30 20];


set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(3,1,[.06 .02],[.08 .05],[.09 .06])

%ha = tight_subplot(3,1,[.02 .00],[.02 .02],[.05 .03])

axes(ha(1))
ax1=gca;

plot_3d_curr_global_quiv(curlx,Xgrid,Ygrid,Zgrid,'\it{J_x}\rm{  (nA m^{-2})}',curlx,curly,curlz)

%plot_3d_ion_global(curlx,Xgrid,Ygrid,Zgrid,'J_x (nA/m^2)')

text(2.05,-1,4.2,'d','Fontsize',20,'FontWeight','bold','FontName','helvetica')

xlabel([])

axes(ha(2))
ax1=gca;

plot_3d_curr_global_quiv(curly,Xgrid,Ygrid,Zgrid,'\it{J_y}\rm{  (nA m^{-2})}',curlx,curly,curlz)

%plot_3d_ion_global(curly,Xgrid,Ygrid,Zgrid,'J_x (nA/m^2)')

text(2.05,-1,4.2,'e','Fontsize',20,'FontWeight','bold','FontName','helvetica')

xlabel([])

axes(ha(3))
ax1=gca;

plot_3d_curr_global_quiv(curlz,Xgrid,Ygrid,Zgrid,'\it{J_z}\rm{  (nA m^{-2})}',curlx,curly,curlz)

%plot_3d_ion_global(curlz,Xgrid,Ygrid,Zgrid,'J_x (nA/m^2)')

text(2.05,-1,4.2,'f','Fontsize',20,'FontWeight','bold','FontName','helvetica')

colormap(mycolormap8)


%%


%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/current/Jmse_xyz_slice_3_0.5_t1107.tiff']);


%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/current/Jmse_xyz_lines_3_0.5_t1105.tiff']);

%print('-dtiff ','-r300',['/Users/gaojiawei/Desktop/global/tu/current/Jmse_xyz_lines_3_0.5_psw-05.tiff']);


print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Jmse_xyz_slice_3_0.5_Bt10_3_q1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Jmse_xyz_slice_3_0.5_Bt10_3_q1.eps']);


