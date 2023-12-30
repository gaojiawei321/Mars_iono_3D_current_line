

draw_data=mola128_dx;

dlevels = [0,10,20,25,50,200,500,800,1000,1100,1200,1500,2000,3000,4000] ;
   
for k = 1 : length(dlevels) - 1
     
   draw_data(find(draw_data>dlevels(k) & draw_data<=dlevels(k+1))) = k ;
   
end

   draw_data(find(draw_data==dlevels(1))) = 1 ;


   cmap = colormap(jet(length(dlevels) - 1)) ;

   
   pcolor(draw_data)
   
   colormap(cmap) ;
    
   caxis([0 length(dlevels)-1]) ;
    
   shading flat
   

   cbar = colorbar ;
    
   set(cbar,'Ticks', 0:1:length(dlevels)-1,'TickLabels',dlevels) ;