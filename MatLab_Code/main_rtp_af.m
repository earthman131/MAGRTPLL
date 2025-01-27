%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The antisymmetric factor method for magnetic reduction to the pole (rtp) at low latitudes
% Author: Lianghui Guo (guolh@cugb.edu.cn)
% Organization: China University of Geosciences (Beijing), School of Geophysics and Information Technology
% Compiled version: MATLAB R2017b
% Reference: 
%       Guo L H, Shi L, Meng X H. The antisymmetric factor method for magnetic reduction
%       to the pole at low latitudes. Journal of Applied Geophysics, 2013, 92, 103-109.
% Description of the input parameters: 
%       infile: original magnetic anomaly data file
%       inc: inclination, unit: ��
%       dec: declination, unit: ��
%       af0: starting Angle of antisymmetric factor method, unit: ��
% Description of the output parameters: 
%       outfile: result file
% Description of primary identifiers��
%       x, y: x, y verctor
%       nx, ny: number of points in x and y directions
%       dx, dy: spacing in x and y directions
%       g: anomaly before rtp
%       npts: extension points
%       grtp: anomaly after rtp
% Description of subroutine function��
%       readgrd.m: read surfer text grd file
%       taper2d.m: cosine attenuation edge extension
%       rtp_af.m: antisymmetric factor method for magnetic reduction to the pole (rtp)
%       savegrd.m: save surfer text grd file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%% I/O parameters %%%%%%%%%%%%%
infile = 'anomaly.grd'; 
outfile = 'anomaly_rtp_af.grd'; 
inc = 0; 
dec = 0; 
af0 = 60; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g,x,y,nx,ny,dx,dy] = readgrd(infile);
[X,Y] = meshgrid(x,y);
nmax = max([nx ny]);
npts = 2^nextpow2(nmax);
% The antisymmetric factor method for magnetic reduction to the pole (rtp)
grtp = rtp_af(g,npts,ny,nx,dec,inc,af0);
% Save
savegrd(grtp,x,y,nx,ny,outfile);