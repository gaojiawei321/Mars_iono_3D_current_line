clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate current in Martian ionosphere
% 2022_10_30 modified JW Gao


%load ionMSEsph_Brtp_L10_-F50_20221113
%load ionMSEsph_Brtp_Bt10_20221109

%load ionMSEsph_Brtp_L10_20221113
%load ionMSEsph_Brtp_L10_Bm150_20221220
load ionMSEsph_Brtp_L10_Bm120_Amse_20230312

%load ionMSOsph_Brtp_L10_Bm120L19_Amse_2023

%load ionMSEsph_Brtp_L10_Bm120_Amse_20221222


%load ionMSEsph_Brtp_L10_Bm150_20221220

load mycolormap9
load myredblue1


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


%%
%save current_MSE_iono_1206 curlx curly curlz xgj ygj zgj

%load /Users/gaojiawei/program/mars/external/suan/Conduc_20130119.mat

%load Conduc_20230207.mat
load Conduc_202309.mat
%load Conduc_2023MSO.mat

%load Conduc_20230705.mat


k = (1/4)*ones(2);

Para1=conv2(Para,k,'same');
Hall1=conv2(Hall,k,'same');
Peder1=conv2(Peder,k,'same');



%% calcualte at J(1,:,:)

for j=1:yge
    for k=1:zge
        
           r=xgj(i);
            theta=ygj(j);
            phi=zgj(k);
            
        ce=1; % altitude ceng 150 km
        normb=[quanBx(ce,j,k),quanBy(ce,j,k),quanBz(ce,j,k)]/norm([quanBx(ce,j,k),quanBy(ce,j,k),quanBz(ce,j,k)]);
        
        J=[curlx(ce,j,k),curly(ce,j,k),curlz(ce,j,k)];
        
        J_para=dot(J,normb)*normb;
        J_perp=J-J_para;
        
        J_para_n(j,k)=norm((J_para));
        J_perp_n(j,k)=norm((J_perp));
        
        E_para(j,k,:)=J_para/Para1(j,k);
        E_perp(j,k,:)=(Peder1(j,k)/(Peder1(j,k)^2+Hall1(j,k)^2)) * (J_perp+ Hall1(j,k)/Peder1(j,k)* cross(J_perp,normb))
        
        E_para_n(j,k)=norm(squeeze(E_para(j,k,:)));
        E_perp_n(j,k)=norm(squeeze(E_perp(j,k,:)));
        
        Ez=E_para(j,k,:)+E_perp(j,k,:);
        E_x(j,k)=Ez(1);
        E_y(j,k)=Ez(2);
        E_z(j,k)=Ez(3);
        
        
  %      E_tot_n(j,k,:)=norm(squeeze(E_tot(j,k,:)));
        
        Ez1=datac2s(E_x(j,k),E_y(j,k),E_z(j,k),pi/2-theta, phi);
        
                    
     %  [Exyz]=datas2c(Ez1(1),Ez1(2),Ez1(3),pi/2-theta, phi);
            
     
        Er(j,k)=Ez1(1);
        Et(j,k)=Ez1(2);
        Ep(j,k)=Ez1(3);
        
    end
end
%%
%Ertp1=Ertp;

%Ertp1(Ertp1>0)=log10(+Ertp1(Ertp1>0));
%Ertp1(Ertp1<0)=-log10(-Ertp1(Ertp1<0));

%E_tot1=E_tot;

%E_tot1(E_tot1>0)=log10(+E_tot1(E_tot1>0));
%E_tot1(E_tot1<0)=-log10(-E_tot1(E_tot1<0));




%% plot electron field


figure
pos=[1 1 40 10];
set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,3,[.03 .03],[.1 .01],[.03 .01])

axes(ha(1))
ax1=gca;



plot_proj_Er_iono(Er,Ygrid,Zgrid,Er,Et,Ep)
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['d'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['E_r'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(2))
ax1=gca;

plot_proj_Er_iono(Et,Ygrid,Zgrid,Er,Et,Ep)
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['e'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['E_\theta'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(3))
ax1=gca;

plot_proj_Er_iono(Ep,Ygrid,Zgrid,Er,Et,Ep)
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['f'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['E_\phi'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

h1=colorbar('southoutside');
h1.Position=[0.3 0.18 0.4 0.02]
h1.Label.String = '\it{E}\rm{  (nV m^{-1})}';
h1.FontSize=12;
h1.FontWeight='bold';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;

h1.TickLabels={'-10^{6}','-10^{4}','-10^{2}','10^{0}','10^{2}','10^{4}','10^{6}'}

%colormap(mycolormap9)
colormap(myredblue1)

%%



print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/F4_Ion_MSE_xyz_E1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/F4_Ion_MSE_xyz_E1.eps']);


%print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F4_Ion_MSE_xyz_E1.tiff']);

%print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F4_Ion_MSE_xyz_E1.eps']);

%%
%% plot cond field


figure
pos=[1 1 40 10];
set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(1,3,[.03 .03],[.1 .01],[.03 .01])

axes(ha(1))
ax1=gca;

plot_proj_Con_iono(log10(Hall1),Ygrid,Zgrid)

%plot_proj_Er_iono(E_tot(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['a'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['Hall'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(2))
ax1=gca;

plot_proj_Con_iono(log10(Peder1),Ygrid,Zgrid)

%plot_proj_Er_iono(E_tot(:,:,2),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['b'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['Pedersen'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(3))
ax1=gca;

plot_proj_Con_iono(log10(Para1),Ygrid,Zgrid)

%plot_proj_Er_iono(E_tot(:,:,3),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))
tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['c'],'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['Parallel'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=150 km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

h1=colorbar('southoutside');
h1.Position=[0.3 0.18 0.4 0.02]
h1.Label.String = '\it{Conductivity}\rm{  (S m^{-1})}';
h1.FontSize=12;
h1.FontWeight='bold';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;

h1.Ticks=[-4,-3,-2,-1,0,1];
h1.TickLabels={'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}'}


%colormap(mycolormap9)
%colormap(myredblue1)
%%

print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/F4_Ion_MSE_xyz_C1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/F4_Ion_MSE_xyz_C1.eps']);


%print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F4_Ion_MSE_xyz_C1.tiff']);

%print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F4_Ion_MSE_xyz_C1.eps']);


%% other

figure
[s1,s2]=size(E_tot(:,:,1));
Eplot=zeros(s1+1,s2+1);
Eplot(1:s1,1:s2)=Ertp(:,:,1);


%axesm mollweid ;
%framem; %gridm
%pcolorm(Ygrid/pi*180,Zgrid/pi*180,(Eplot)) ;  hold on;
%pcolor(log(Eplot))

plot_proj_Er_iono(E_tot(:,:,1),Ygrid,Zgrid,Ertp(:,:,1),Ertp(:,:,2),Ertp(:,:,3))
%plot_proj_Er_iono(Erpt(:,:,1),Ygrid,Zgrid,Erpt(:,:,1),Erpt(:,:,2),Erpt(:,:,3))


title('log ( E_tot_x nV/m) ' )

%cun1=['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/ion/E/E_tot_x_20230130.tiff']
%print('-dtiff ','-r300',cun1);
colorbar

colormap(mycolormap9)




figure
pcolor(log(E_perp_n))

title('log ( E_perp_n ) nV/m' )

%cun1=['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/ion/E/E_perp_n_20230119.tiff']
print('-dtiff ','-r300',cun1);
colorbar

a=1




