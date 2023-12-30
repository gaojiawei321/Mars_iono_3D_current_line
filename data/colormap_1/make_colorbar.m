clear all; close all;clc

load mycolormap7
load mycolormap8
load mycolormap9

load mycolormap10
load mycolormap1


colorbar
colormap(mycolormap1)

cmap1=get(gca,'colormap')

mycolormapd1=cmap1;

save mycolormapd1 mycolormapd1


%save mycolormap8 mycolormap8

