clear;clc;close all;
data=xlsread('log4D.csv');
%data= importdata('log4D.dat');
%print data;
x=data(:,2);y=data(:,1);z=data(:,3);e=data(:,4);
%figure('Renderer','zbuffer','Color',[1 1 1]);
deltax=linspace(min(x),max(x),50);
deltay=linspace(min(y),max(y),50);
deltaz=linspace(min(z),max(z),20);
[xq,yq,zq] = meshgrid(deltax, deltay,deltaz);
vq = griddata(x,y,z,e,xq,yq,zq,'natural');

xs = deltax;
ys = deltay;
zs = deltaz;

h = slice(xq,yq,zq,vq,xs,ys,zs);
set(h,'FaceColor','interp',...
    'EdgeColor','none')
camproj perspective
box on
colormap jet
colorbar
