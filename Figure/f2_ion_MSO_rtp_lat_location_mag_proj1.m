clear all;close all;clc

%}

load ionMSOsph_Brtp_L10_Bm120_2023
%load ionMSOsph_Brtp_Un_2023

load mycolormap9
load myredblue
yge=10;
zge=20;

% PCOLOR
%qz=permute(quanBp,[2 1 3]);


%% Tu iono
%% Tu iono

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
        tp1=['plot_proj_mag_iono(quanBr(ceng',num2str(hc),',:,:),Ygrid,Zgrid,quanBr(ceng',num2str(hc),',:,:),quanBt(ceng',num2str(hc),',:,:),quanBp(ceng',num2str(hc),',:,:))'];
        
        eval(tp1)
    else if hcj==2
            tp1=['plot_proj_mag_iono(quanBt(ceng',num2str(hc),',:,:),Ygrid,Zgrid,quanBr(ceng',num2str(hc),',:,:),quanBt(ceng',num2str(hc),',:,:),quanBp(ceng',num2str(hc),',:,:))'];
            
            eval(tp1)
        else
            tp1=['plot_proj_mag_iono(quanBp(ceng',num2str(hc),',:,:),Ygrid,Zgrid,quanBr(ceng',num2str(hc),',:,:),quanBt(ceng',num2str(hc),',:,:),quanBp(ceng',num2str(hc),',:,:))'];
            
            eval(tp1)
        end
    end
    
    
    tj1=-2; tj2=1.55; tj3=1.38;    tj4=1.52;
    
    list1=['abcdefghijklmnopqrstuvwxyz'];
    
    
    text(-3.2,1.8,list1(i),'FontSize',20,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;
    
    
    if hcj==1
        
        text(tj1,tj2,['B_r'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
    else if hcj==2
            text(tj1,tj2,['B_\theta'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
            
        else
            text(tj1,tj2,['B_\phi'],'FontSize',15,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;
            
        end
    end
    
    t2=['ceng',num2str(hc)];
    
    text(tj3,tj4,['h=',num2str(Ch( eval(t2)  )),' km'],'FontSize',12,'fontweight', 'bold','FontName','times'); hold on;
    
    
end

h1=colorbar('southoutside');
h1.Position=[0.3 0.18 0.4 0.02]

h1.Label.String = '\it{B}\rm{  (nT)}';
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

plot_proj_mag_iono(quanBr(ceng1,:,:),Ygrid,Zgrid,quanBr(ceng1,:,:),quanBt(ceng1,:,:),quanBp(ceng1,:,:))

% title('J_x','FontSize',16)

%VerticalAlignment

tj1=-1.8; tj2=1.56; tj3=1.38;    tj4=1.52;

d1=25;
d2=20;

text(-3.4,1.8,['a'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_r'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(2))

plot_proj_mag_iono(quanBt(ceng1,:,:),Ygrid,Zgrid,quanBr(ceng1,:,:),quanBt(ceng1,:,:),quanBp(ceng1,:,:))


text(-3.4,1.8,['b'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_\theta'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(3))

plot_proj_mag_iono(quanBp(ceng1,:,:),Ygrid,Zgrid,quanBr(ceng1,:,:),quanBt(ceng1,:,:),quanBp(ceng1,:,:))


text(-3.4,1.8,['c'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_\phi'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng1)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


axes(ha(4))

plot_proj_mag_iono(quanBr(ceng2,:,:),Ygrid,Zgrid,quanBr(ceng2,:,:),quanBt(ceng2,:,:),quanBp(ceng2,:,:))

text(-3.4,1.8,['d'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_r'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(5))

plot_proj_mag_iono(quanBt(ceng2,:,:),Ygrid,Zgrid,quanBr(ceng2,:,:),quanBt(ceng2,:,:),quanBp(ceng2,:,:))


text(-3.4,1.8,['e'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_\theta'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;

axes(ha(6))
ax1=gca;

plot_proj_mag_iono(quanBp(ceng2,:,:),Ygrid,Zgrid,quanBr(ceng2,:,:),quanBt(ceng2,:,:),quanBp(ceng2,:,:))

text(-3.4,1.8,['f'],'FontSize',d1,'fontweight', 'bold','FontName','helvetica','HorizontalAlignment','center'); hold on;

text(tj1,tj2,['B_\phi'],'FontSize',d2,'FontAngle','italic','fontweight', 'bold','FontName','times'); hold on;

text(tj3,tj4,['h=',num2str(Ch(ceng2)),' km'],'FontSize',15,'fontweight', 'bold','FontName','times'); hold on;


h1=colorbar;
h1.Position=[0.925 0.1 0.015 0.8]
h1.Label.String = '\it{B}\rm{  (nT)}';
h1.FontSize=14;
h1.FontWeight='normal';
h1.TickDirection='out';
h1.TickLength=0.01;
h1.LineWidth=1;

load mycolormap12
colormap(mycolormap12)
%colormap(myredblue)

%text(-1.7,-2.0,'10 nT','Rotation',0,'FontSize',20); hold on;

%text(-1.3,-2.3,'\rightarrow','Rotation',0,'FontSize',20); hold on;



%%
print('-dtiff','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F2_Ion_MSO_rtp_Bm10_9_q1.tiff']);

print('-depsc','-r300',['/Users/gaojiawei/Desktop/Baiducloud/work/current/tu/Paper_tu/Fig2/F2_Ion_MSO_rtp_Bm10_9_q1.eps']);


print('-dtiff','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F2_Ion_MSE_rtp_Bm10_9_q1.tiff']);

print('-depsc','-r300',['C:\Users\gao\Desktop\BaiduSyncdisk\work\current\tu\Paper_tu\F2_Ion_MSE_rtp_Bm10_9_q1.eps']);


a=1

%%