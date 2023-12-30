clear all;close;clc

load mycolormap1
load mycolormap2
load mycolormap3
load mycolormap4

load mycolormap6

load mycolormap7
load mycolormap8
load mycolormap9

load mycolormap10

a=[1 1;1 1]
pcolor(a)

colormap(mycolormap10)

%colormap(hot)
colorbar


cmap1=get(gca,'colormap')

%cmap2=flip(cmap1)
%colormap(cmap2)

mycolormap11=cmap1;


save mycolormap11 mycolormap11
