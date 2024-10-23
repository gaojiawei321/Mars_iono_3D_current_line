clear all;close all;clc

%% model at 150 km


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


%%



%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate current in Martian ionosphere
% 2022_11_13

%
%load MAVEN_4s_q_2022.mat

load /Users/gaojiawei/program/mars/external/quan/MAVEN_4s_f_2022.mat


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


%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 );

hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmso(:,1)~=0 & Bmodel150<10 );

%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 & MAVEN_4s_f.mse_pan==1);

%hou=find( L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0  & Bmodel150<10 & MAVEN_4s_f.mse_pan==1 & MAVEN_4s_f.F107>50);

% add z hemisphere
%hou=find(Bmt<10 & L>l1& L<l2 & MAVEN_4s_f.Lmse_zh(:,1)~=0 & MAVEN_4s_f.Lpc(:,3)>0 );

%[La,Le,Lr] = cart2sph( MAVEN_4s_f.Lmse_zh(hou,1),MAVEN_4s_f.Lmse_zh(hou,2),MAVEN_4s_f.Lmse_zh(hou,3));


[La,Le,Lr] = cart2sph( MAVEN_4s_f.Lmso(hou,1),MAVEN_4s_f.Lmso(hou,2),MAVEN_4s_f.Lmso(hou,3));


Let=pi/2-Le;

%B in MSE; cart to sph
%% check this on: seem make sence


%check complete



Bmse_sph=datac2s(MAVEN_4s_f.Bmso(hou,1),MAVEN_4s_f.Bmso(hou,2),MAVEN_4s_f.Bmso(hou,3),Let',La');


a=1

%%%%%% Bmse sph (BR BT BP)   Location mse_sph (LR LT LP)
%%%%%%  Br(ourward) Btheta(southward) Bphi(eastward)

bsph=[ Bmse_sph', Lr, Le, La];

bxyz=[ MAVEN_4s_f.Bmso(hou,:), Lr, Le, La];

Ntene=[ MAVEN_4s_f.KPnco2(hou,:),MAVEN_4s_f.te(hou,:),MAVEN_4s_f.ne(hou,:), Lr, Le, La];
% Nco2 Te Ne
z1=find(Ntene==0);
Ntene(z1)=NaN;

%clear MAVEN_4s_f


%%%%
%z1=find(b(:,1)==0);
%b(z1,:)=[];

xge=10;
yge=10;
zge=20;

xs=125+Rm;
xe=625+Rm;

ys=-pi/2;
ye=pi/2;

zs=-pi;
ze=pi;

xj=(xe-xs)/xge;
yj=(ye-ys)/yge;
zj=(ze-zs)/zge;


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

quanN1=zeros(xge,yge,zge);
quanN2=zeros(xge,yge,zge);
quanN3=zeros(xge,yge,zge);

ge=zeros(xge,yge,zge);

for i=1:xge
    for j=1:yge
        for k=1:zge

            i
            j
            k


            z1=find(bsph(:,4)<Xgrid(i+1) &bsph(:,4)>Xgrid(i)    &bsph(:,5)<Ygrid(j+1)&bsph(:,5)> Ygrid(j)  &bsph(:,6)<Zgrid(k+1)  &bsph(:,6)>Zgrid(k) ) ;

            %median

       
            Bmse_sph=datac2s(MAVEN_4s_f.Bmso(hou,1),MAVEN_4s_f.Bmso(hou,2),MAVEN_4s_f.Bmso(hou,3),Let',La');

       
            quanBr(i,j,k)=nanmean(bsph(z1,1));
            quanBt(i,j,k)=nanmean(bsph(z1,2));
            quanBp(i,j,k)=nanmean(bsph(z1,3));

            quanBx(i,j,k)=nanmean(bxyz(z1,1));
            quanBy(i,j,k)=nanmean(bxyz(z1,2));
            quanBz(i,j,k)=nanmean(bxyz(z1,3));

   
            quanN1(i,j,k)=nanmean(Ntene(z1,1));
            quanN2(i,j,k)=nanmean(Ntene(z1,2));
            quanN3(i,j,k)=nanmean(Ntene(z1,3));

            fn1=find(Ntene(z1,3)<1);
            Ntene(z1(fn1),3)=NaN;
            fn1=find(Ntene(z1,2)<1);
            Ntene(z1(fn1),2)=NaN;

            Q= quantile(Ntene(z1,1),[0.25,0.5, 0.75]) ;
            QN1l(i,j,k)=Q(1);
            QuanN1(i,j,k)=Q(2);
            QN1u(i,j,k)=Q(2);
            Q= quantile(Ntene(z1,2),[0.25, 0.5, 0.75]) ;
            QN2l(i,j,k)=Q(1);
            QuanN2(i,j,k)=Q(2);
            QN2u(i,j,k)=Q(2);
            Q= quantile(Ntene(z1,3),[0.25,0.5, 0.75]) ;
            QN3l(i,j,k)=Q(1);
            QuanN3(i,j,k)=Q(2);
            QN3u(i,j,k)=Q(3);

            %  ttest

            %       [h,p] = ttest(bxyz(z1,1));
            %       [h,p] = ttest(bxyz(z1,2));
            %       [h,p,ci,stats]  = ttest(bxyz(z1,3));


            ge_sph(i,j,k)=length(bsph(z1,3));
            ge_xyz(i,j,k)=length(bxyz(z1,3));
            ge_N3(i,j,k)=length(Ntene(z1,3));
            
            
            
            
            
        end

    end
end


  for i=1:xge
    for j=1:yge
        for k=1:zge     
                  [lx,ly,lz] = sph2cart(Zgrid(k),Ygrid(j),Xgrid(i));

           
            r=(lx^2+ly^2+lz^2)^0.5;
            SZA=acos(lx/r);
            Nmax(i,j,k)=2.58 * 10^4 * (50* cos(SZA))^0.5;
           end

    end
end
        

%%
save Ndata2023_MSO quanBp quanBr quanBt quanN1 quanN2 quanN3 QuanN1 QuanN2 QuanN3  QN1l QN1u QN2l QN2u QN3l QN3u

%check N Te Ne data
%s1=squeeze(nanmean(quanN1,1));
squeeze(ge_sph(1,:,:));

%real(squeeze(Nmax(1,:,:)) );

%pcolor(s1)

method='linear';
quanBx=fillmissing(quanBx,method);
quanBy=fillmissing(quanBy,method);
quanBz=fillmissing(quanBz,method);


quanBr=fillmissing(quanBr,method);
quanBt=fillmissing(quanBt,method);
quanBp=fillmissing(quanBp,method);

QuanN2=fillmissing(QuanN2,method);

QuanN3=fillmissing(QuanN3,method);

s1=squeeze(QuanN1(1,:,:));
%pcolor(s1)
%caxis([0 1*10^10])

Te1=squeeze(QuanN2(1,:,:));
Ne1=squeeze(QuanN3(1,:,:));



TeMSE=load('Ndata2023.mat');
TeE=squeeze(TeMSE.QuanN2(1,:,:));

NeE=squeeze(TeMSE.QuanN3(1,:,:));


%% 150 km Mars

e=1.6*10^-19;
me=9.1*10^-31;
mi=1.67*10^-27; 


%% 剖面

%h=100:40:400;
%h=200;

%Ne=[50000 100000 50000 30000 20000 15000 10000  5000] *10^6;
Ne=Ne1 *10^6;

Ne=real(squeeze(Nmax(1,:,:)) ) *10^6;

Ne=NeE *10^6;


%Te=[100 120 400 1600 2000 2700 3000 3300]
%Ti=[100 120 160 400 1000 2700 3000 3300]

%Tr=(Te+Ti)/2

%Nco2=[1*10^12 1*10^10 1*10^8 1*10^6 1*10^5 5*10^4 3*10^4 1*10^4] *10^6;
Nco2=[1*10^10 ] *10^6;

%Nco2=[nanmean(nanmean(s1)) ] *10^6;

%No

%Ven=3.68*10^-8* Nco2/10^6.* (1+4.1*10^-11.* (abs(4500-Te)).^2.93 )

%s1
Ven1=3.68*10^-8* Nco2/10^6.* (1+4.1*10^-11.* (abs(4500-TeE)).^2.93 )


% pcolor(Ven1)
% caxis([500 5000])

%用氧的。 O+,O
%Vi=3.67*10^-11*Nco2./10^6.*Tr.^0.5.*(1-0.064*log10(Tr)).^2


%用氧的。 H+,O
%Vi=6.61*10^-11*Nco2./10^6.*Ti.^0.5.*(1-0.047*log10(Ti)).^2

% nonresonant collision frequency。O+2 CO2

%Vin=(5.63*10^-10)*Nco2*10^-6;

%s1
Vin=(5.63*10^-10)*Nco2*10^-6;

Vin1=Vin*ones(size(Ne1));

%Vin=Vin+Vi;



%% Vin（Nco2）
%% Ven（Nco2, Te）

%% Pare: Ne Ven B
%% Peder: Ne Vin B

%% Nco2 Te Ne B


%Vei=[1000000 200 10   1     1*10^-2 1*10^-3    1*10^-5  1*10^-6];


%B=[6 6 30 50 50 50 50 50] *10^-9;




B=(quanBr.^2+quanBt.^2+quanBp.^2).^0.5;
B1=squeeze(B(1,:,:))*10^-9;

%B=[5 5 5 5 5 5 5 5] *10^-9;

Fe=e*B1/me;
Fi=e*B1/(mi*16);


%Para1=Ne*e^2/me./Ven

%Peder1=Ne*e^2.*Ven/me./( Ven.^2+  Fe.^2   )  ;

%Hall1=Ne*e^2.*Fe/me./( Ven.^2+  Fe.^2   )  ;


qmi=32;

Para=Ne*e^2/me./(Ven1);


%Peder=Ne*e^2.* (Vi/(mi*qmi)./( Vi.^2+  Fi.^2 )   +  Ven/me./( Ven.^2+  Fe.^2 ) ) ;

%Hall=Ne*e^2.*  ( Fe/me./( Ven.^2+  Fe.^2 ) -  Fi/(mi*qmi)./( Vi.^2+  Fi.^2 ) ) ;

Peder=Ne*e^2.* (Vin1/(mi*qmi)./( Vin1.^2+  Fi.^2 )   +  Ven1/me./( Ven1.^2+  Fe.^2 ) ) ;

Hall=Ne*e^2.*  ( Fe/me./( Ven1.^2+  Fe.^2 ) -  Fi/(mi*qmi)./( Vin1.^2+  Fi.^2 ) ) ;


%save Conduc_20230207 Para Peder Hall
save Conduc_2023MSO Para Peder Hall

%%







%%