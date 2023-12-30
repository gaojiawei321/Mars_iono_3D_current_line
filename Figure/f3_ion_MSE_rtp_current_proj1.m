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

load ionMSEsph_Brtp_L10_Bm120_Amse_2023

%load ionMSOsph_Brtp_L10_Bm120L19_Amse_2023

%load ionMSEsph_Brtp_L10_Bm120_Amse_20221222


%load ionMSEsph_Brtp_L10_Bm150_20221220

load mycolormap9


%%
%method='linear';
%quanBx=fillmissing(quanBx,method);
%quanBy=fillmissing(quanBy,method);
%quanBz=fillmissing(quanBz,method);


%quanBr=fillmissing(quanBr,method);
%quanBt=fillmissing(quanBt,method);
%quanBp=fillmissing(quanBp,method);


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
            
       %      JJ=datac2s(Jxyz(1),Jxyz(2),Jxyz(3),pi/2-theta, phi);
            
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

%%


%% Tu iono proj9
%{


figure
pos=[1 1 40 22];


set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(3,3,[.03 .03],[.12 .04],[.03 .03])

ceng1=7;
ceng2=4;
ceng3=1;
Ch=((Xgrid(2:end)-3393.5)+(Xgrid(1:end-1)-3393.5))/2;

axes(ha(1))
ax1=gca;

plot_proj_curr_iono(curlx(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))

% title('J_x','FontSize',16)

%VerticalAlignment
d1=25;
d2=20;

tj1=-1.8; tj2=1.52; tj3=1.38;    tj4=1.52;

text(-3.4,1.8,['j'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_x'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(2))


plot_proj_curr_iono(curly(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))


text(-3.4,1.8,['k'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_y'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(3))

plot_proj_curr_iono(curlz(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))

text(-3.4,1.8,['l'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_z'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(4))

plot_proj_curr_iono(curlx(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))


text(-3.4,1.8,['m'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_x'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(5))

plot_proj_curr_iono(curly(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))


text(-3.4,1.8,['n'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_y'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(6))
ax1=gca;

plot_proj_curr_iono(curlz(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))


text(-3.4,1.8,['o'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_z'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(7))

plot_proj_curr_iono(curlx(ceng3,:,:),Ygrid,Zgrid, curlr(ceng3,:,:),curlt(ceng3,:,:),curlp(ceng3,:,:))

text(-3.4,1.8,['p'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_x'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng3)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(8))

plot_proj_curr_iono(curly(ceng3,:,:),Ygrid,Zgrid, curlr(ceng3,:,:),curlt(ceng3,:,:),curlp(ceng3,:,:))

text(-3.4,1.8,['q'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_y'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng3)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(9))

plot_proj_curr_iono(curlz(ceng3,:,:),Ygrid,Zgrid, curlr(ceng3,:,:),curlt(ceng3,:,:),curlp(ceng3,:,:))


text(-3.4,1.8,['r'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_z'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng3)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


h1=colorbar('southoutside');
h1.Position=[0.3 0.09 0.4 0.02]
h1.Label.String = '\it{J}\rm{  (nA m^{-2})}';
h1.FontSize=15;
h1.FontWeight='normal';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;

colormap(mycolormap9)


%text(-1.7,-2.0,'50 nA m^{-2}','Rotation',0,'FontSize',20); hold on;

%text(-1.3,-2.3,'\rightarrow','Rotation',0,'FontSize',20); hold on;

cun1=['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/ion/J/Ion_MSE_xyz_curr_Bm10_9_q1.tiff']

%}

%%


figure
pos=[1 1 40 15];


set(gcf,'unit','centimeters','position',pos)

ha = tight_subplot(2,3,[.03 .03],[.04 .04],[.03 .09])

ceng1=3;
ceng2=1;
%ceng3=1;
Ch=((Xgrid(2:end)-3393.5)+(Xgrid(1:end-1)-3393.5))/2;

axes(ha(1))
ax1=gca;

plot_proj_curr_iono(curlr(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))

% title('J_x','FontSize',16)

%VerticalAlignment

tj1=-1.8; tj2=1.56; tj3=1.38;    tj4=1.52;

d1=25;
d2=20;

text(-3.4,1.8,['g'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_r'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(2))

plot_proj_curr_iono(curlt(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))


text(-3.4,1.8,['h'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_\theta'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(3))

plot_proj_curr_iono(curlp(ceng1,:,:),Ygrid,Zgrid, curlr(ceng1,:,:),curlt(ceng1,:,:),curlp(ceng1,:,:))

text(-3.4,1.8,['i'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_\phi'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(4))

plot_proj_curr_iono(curlr(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))

text(-3.4,1.8,['j'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_r'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(5))

plot_proj_curr_iono(curlt(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))


text(-3.4,1.8,['k'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_\theta'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(6))
ax1=gca;

plot_proj_curr_iono(curlp(ceng2,:,:),Ygrid,Zgrid, curlr(ceng2,:,:),curlt(ceng2,:,:),curlp(ceng2,:,:))

text(-3.4,1.8,['l'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['J_\phi'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


h1=colorbar;
h1.Position=[0.925 0.1 0.015 0.8]
h1.Label.String = '\it{J}\rm{  (nA m^{-2})}';
h1.FontSize=14;
h1.FontWeight='normal';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;


%colormap(mycolormap9)
%colormap(myredblue)

%text(-1.7,-2.0,'10 nT','Rotation',0,'FontSize',20); hold on;

%text(-1.3,-2.3,'\rightarrow','Rotation',0,'FontSize',20); hold on;

load mycolormap9
colormap(mycolormap9)

%%


print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F3_Ion_MSE_xyz_curr_Bm10_9_q1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F3_Ion_MSE_xyz_curr_Bm10_9_q1.eps']);




print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F3_Ion_MSE_xyz_curr_Bm10_9_q1.tiff']);

print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F3_Ion_MSE_xyz_curr_Bm10_9_q1.eps']);




a=1
