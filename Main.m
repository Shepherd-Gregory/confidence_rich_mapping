%{  
    Copyright (C) 2021  Yang Xu (xuyang94@zju.edu.cn)
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details. 
%}
clear all;
close all;
%% Params setup
param.mapSize = 33; % in meter
% 1. Decide map resolution, i.e., the number of grids for 1 meter.
param.resol = 3; 
% 2. Decide the initial map size in pixels
param.size = [param.mapSize * param.resol, param.mapSize * param.resol];
% 3. Indicate where you will put the origin in pixels
param.origin = [0,0]'; 
% 4. Log-odd parameters 
param.lo_occ = 1;
param.lo_free = 0.5; 
param.lo_max = 100;
param.lo_min = -100;
param.res_z = 10/param.resol ;
% params for mixture sensor model 
param.z_hit = 0.7;
param.z_max = 0.1;
param.z_rand = 0.1;
param.z_short = 0.1;
param.lambda_short = 0.2;
% integration resolution for map occupancy
param.res_m = 10;
% sensor noise
param.dev = 0.1 / param.resol; 
% Max sensing range
param.maxrange = 8; % in meter
param.flag = 0; % 0 for ijrr map; 1 for wall
%% load your dataset
% load datasets\wallData.mat
load datasets\ijrrCRM.mat

test_scan = [1:60];
test_step = 1;
test_step_end = 1144; 
ts = tic;
[myMap, com, bin_m, bel_m, exp_m] = CRM(scanRanges(test_scan,test_step:test_step_end), scanAngles(test_scan), pose(:,1+test_step:test_step_end+1), param);

t = toc(ts);

% The final grid map: 
figure
imagesc(myMap); 
hold on
colormap('gray'); axis equal;
title('Log-odds')

figure
imagesc(bin_m); grid on
colormap('gray'); axis equal;
title('Estimated binary occupancy')

figure
imagesc(exp_m);
colormap('gray'); axis equal;
title('Mean estimated occupancy')

figure
imagesc(com);
colormap('gray'); axis equal;
title('occupancy from Log-odds')