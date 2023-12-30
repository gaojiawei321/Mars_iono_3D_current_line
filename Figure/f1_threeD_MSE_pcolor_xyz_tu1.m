clear all;close all;clc


%load MAVEN_4s_q.mat

load MAVEN_4s_q_2023.mat

%load MAVEN_4s_f.mat

Rm=3393.5;


%[theta,rho,z] = cart2pol(MAVEN_4s_f.Lmse_zh(:,3),MAVEN_4s_f.Lmse_zh(:,2),MAVEN_4s_f.Lmse_zh(:,1));

%bpol=datac2p(MAVEN_4s_f.Bmse_zh(:,3),MAVEN_4s_f.Bmse_zh(:,2),MAVEN_4s_f.Bmse_zh(:,1),theta')';


%b=[ bpol,MAVEN_4s_f.Lmse_zh./Rm];



%Bimf=vecnorm(MAVEN_4s_f.IMFB')';
%Vimf=vecnorm(MAVEN_4s_f.IMFV')';

Rm=3393.5;
L=( MAVEN_4s_f.Lmso(:,1).^2+MAVEN_4s_f.Lmso(:,2).^2+MAVEN_4s_f.Lmso(:,3).^2).^0.5-Rm;

Bmodel=MAVEN_4s_f.Bmso-MAVEN_4s_f.Bmso_res_g;

Bmt=(Bmodel(:,1).^2+Bmodel(:,2).^2+Bmodel(:,3).^2).^0.5;
Bot=vecnorm(MAVEN_4s_f.Bmso',2);

l1=0;
l2=10000;
%hou=find(Bmt<10);


Mi=1.67*10^-27;
u0=4*pi*10^-7;

%Vsw=vecnorm(MAVEN_4s_f.IMFV',2);
%Psw=   (MAVEN_4s_f.Dsw *10^6 ).*Mi.*(Vsw'*1000).^2 *10^9;
%Bt=vecnorm(MAVEN_4s_f.IMFB',2);
%cone=acos( - (MAVEN_4s_f.IMFB(:,1))./(Bt') );


b=[ MAVEN_4s_f.Bmse_zh(:,1:3),MAVEN_4s_f.Lmse_zh./Rm];
%hou=find(Bmt<10 & L>l1& L<l2);

% MSE steady and Bmodel<10 nt
hou=find(Bmt<10  & MAVEN_4s_f.mse_pan==1);
% MSE steady and 10*Bmodel<Bobs 
%hou=find(10.*Bmt<Bot'  & MAVEN_4s_f.mse_pan==1);


%solar wind condition
%hou=find(Bmt<10 & L>l1& L<l2  & Vimf>400);
%hou=find(Bmt<10 & L>l1& L<l2  & Psw<0.5);
%hou=find(Bmt<10 & L>l1& L<l2  & cone>120/180*pi);



%hou=find(Bmt<10 & L>l1& L<l2  & Vimf>400);

%hou=find(Bmt<10 & L>l1& L<l2 &  Vimf<400& Vimf>0.1);

%hou=find(Bmt<10 & L>l1& L<l2  & Bimf<2& Bimf>0.1);

%hou=find(Bmt<10 & L>l1& L<l2  & Bimf>2);

b=b(hou,:);



 %%
z1=find(b(:,1)==0);
b(z1,:)=[];

xge=30;
yge=30;
zge=30;

xs=-3;
xe=3;

ys=-3;
ye=3;

zs=-3;
ze=3;

xj=(xe-xs)/xge;
yj=(ye-ys)/yge;
zj=(ze-zs)/zge;

Xgrid=linspace(xs,xe,xge+1);
Ygrid=linspace(ys,ye,yge+1);
Zgrid=linspace(zs,ze,zge+1);

xgj=Xgrid(1:end-1)+xj/2;
ygj=Ygrid(1:end-1)+yj/2;
zgj=Zgrid(1:end-1)+zj/2;

quanBx=zeros(xge,yge,zge);
quanBy=zeros(xge,yge,zge);
quanBz=zeros(xge,yge,zge);

ge=zeros(xge,yge,zge);

%%
%{
xge=20;
yge=20;
zge=20;

xs=-2;
xe=2;

ys=-2;
ye=2;

zs=-2;
ze=2;

xj=(xe-xs)/xge;
yj=(ye-ys)/yge;
zj=(ze-zs)/zge;

% -1.9: 0.2 : 2.1

% -2: 0.2 : 2
%}
for i=1:xge
    for j=1:yge
        for k=1:zge
        
        i
        j
        k
        
   %     xf=i/ (xge/(xe-xs))- xe +  (xe-xs)/xge/2 ;
   %     yf=j/ (yge/(ye-ys))- ye +  (ye-ys)/yge/2 ;
   %     zf=k/ (zge/(ze-zs))- ze +  (ze-zs)/zge/2 ;
        
      z1=find(b(:,4)<Xgrid(i+1) &b(:,4)>Xgrid(i)    &b(:,5)<Ygrid(j+1)&b(:,5)> Ygrid(j)  &b(:,6)<Zgrid(k+1)  &b(:,6)>Zgrid(k) ) ;
         
        
    %   z1=find(b(:,4)<xf&b(:,4)>(xf- (xe-xs)/xge   )     &b(:,5)<yf&b(:,5)>(yf- (ye-ys)/yge)   &b(:,6)<zf&b(:,6)>(zf- (ze-zs)/zge) ) ;
         
        %median
        
        quanx(i,j,k,:)=nanmean(b(z1,1));
        quany(i,j,k,:)=nanmean(b(z1,2));
        quanz(i,j,k,:)=nanmean(b(z1,3));
        
        ge(i,j,k,:)=length(b(z1,3));
        end

    end
end


%%


%Xj=[-1.8:0.2:2];

readme='MSE quanx Bx   quanx By  quanz Bz  '
save threeDMSE_B_xyz_Bt10_2023 quanx quany quanz readme ge Xgrid Ygrid Zgrid

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

 
