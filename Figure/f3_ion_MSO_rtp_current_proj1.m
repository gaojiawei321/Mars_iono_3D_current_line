clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%

%load ionMSOsph_Brtp_Un_2023
load ionMSOsph_Brtp_L10_Bm120_2023


load mycolormap9


%%
Rm=3393.5;

xgj=Xgrid(1:end-1)+((Xgrid(2)-Xgrid(1)))/2;  % r grid
xgj=xgj*1000; % km --> m
ygj=Ygrid(1:end-1)+((Ygrid(2)-Ygrid(1)))/2;  % theta grid
zgj=Zgrid(1:end-1)+((Zgrid(2)-Zgrid(1)))/2;  % phi grid

[curlr,curlt,curlp] = ion_cal_curr_sph(quanBr,quanBt,quanBp,xgj,ygj,zgj);


%% Tu iono proj9

figure

pos=[1 1 35 8];

set(gcf,'unit','centimeters','position',pos)

%ha = tight_subplot(1,3,[-.01 .00],[.05 .00],[.00 .00])
ha = tight_subplot(1,3,[.03 .03],[.1 .01],[.03 .01])


ceng1=1;



Ch=((Xgrid(2:end)-3393.5)+(Xgrid(1:end-1)-3393.5))/2;

%%
for i=1:3
    
    axes(ha(i))
    
    ax1=gca;
    
    hc=floor(i/3+0.9);
    hcj=i-(floor(i/3+0.9)-1)*3;
    
  if hcj==1
        tp1=['plot_proj_curr_iono(curlr(ceng',num2str(hc),',:,:),Ygrid,Zgrid,curlr(ceng',num2str(hc),',:,:),curlt(ceng',num2str(hc),',:,:),curlp(ceng',num2str(hc),',:,:))'];
        
        eval(tp1)
    else if hcj==2
            tp1=['plot_proj_curr_iono(curlt(ceng',num2str(hc),',:,:),Ygrid,Zgrid,curlr(ceng',num2str(hc),',:,:),curlt(ceng',num2str(hc),',:,:),curlp(ceng',num2str(hc),',:,:))'];
            
            eval(tp1)
        else
            tp1=['plot_proj_curr_iono(curlp(ceng',num2str(hc),',:,:),Ygrid,Zgrid,curlr(ceng',num2str(hc),',:,:),curlt(ceng',num2str(hc),',:,:),curlp(ceng',num2str(hc),',:,:))'];
            
            eval(tp1)
        end
    end
    
    
    tj1=-2; tj2=1.55; tj3=1.38;    tj4=1.52;
    
    %list1=['abcdefghijklmnopqrstuvwxyz'];
    
    list1=['defghijklmnopqrstuvwxyz'];
    
    
    text(-3.2,1.8,list1(i),'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;
    
    
  if hcj==1
        
        text(tj1,tj2,['J_r'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
    else if hcj==2
            text(tj1,tj2,['J_\theta'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
            
        else
            text(tj1,tj2,['J_\phi'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
            
        end
    end
    
  
    
    t2=['ceng',num2str(hc)];
    
    text(tj3,tj4,['h=',num2str(Ch( eval(t2)  )),' km'],'FontSize',12,'fontweight', 'bold','FontName','times'); hold on;
    
    
end

h1=colorbar('southoutside');
h1.Position=[0.3 0.18 0.4 0.02]
h1.Label.String = '\it{J}\rm{  (nA m^{-2})}';
h1.FontSize=12;
h1.FontWeight='bold';

h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;

colormap(mycolormap9)

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


print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F3_Ion_MSO_xyz_curr_Bm10_9_q1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F3_Ion_MSO_xyz_curr_Bm10_9_q1.eps']);




print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F3_Ion_MSO_xyz_curr_Bm10_9_q1.tiff']);

print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F3_Ion_MSO_xyz_curr_Bm10_9_q1.eps']);



a=1
