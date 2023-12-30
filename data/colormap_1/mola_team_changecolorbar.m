clear;close;clc


load mola128_oc
load molaTeam

[a1,a2]=size(mola128_oc)
x=1:1:a1;
y=1:1:a2;
xx=mola128_oc(x,y);
dx=double(xx);
dx=flip(dx);
%mola128_dx=dx;

draw_data=dx;
draw_data1=dx;


dlevels = molaTeam(:,1);

for k = 1 : length(dlevels) - 1
     
    z1=find(draw_data>dlevels(k) & draw_data<=dlevels(k+1));
   draw_data1(z1) = k ;
   
end

   draw_data1(find(draw_data==dlevels(1))) = 1 ;

   
   cmap = molaTeam(:,2:4)/255;

 %  cmap = colormap(jet(length(dlevels) - 1)) ;

   figure
   
   pcolor(draw_data1)
   
   colormap(cmap) ;
    
   caxis([0 length(dlevels)-1]) ;
    
   shading flat
   

   cbar = colorbar ;
    
   set(cbar,'Ticks', 0:1:length(dlevels)-1,'TickLabels',dlevels) ;