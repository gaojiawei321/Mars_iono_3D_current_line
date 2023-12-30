function [curlr,curlt,curlp] = ion_cal_curr_sph(quanBr,quanBt,quanBp,xgj,ygj,zgj)
%ION_CAL_CURR_SPH 此处显示有关此函数的摘要
%   此处显示详细说明

Br=quanBr.*10^-9;
Bt=quanBt.*10^-9;
Bp=quanBp.*10^-9;

[xge,yge,zge]=size(quanBr);


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


end

