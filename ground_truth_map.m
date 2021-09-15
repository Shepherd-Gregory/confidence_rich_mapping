function Ground_truth_map = ground_truth_map(res, mapsize, FLAG)
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
r = res;
size = mapsize;

Ground_truth_map = zeros(size * r, size * r);
if FLAG == 0 
    %% IJRR MAP
    Ground_truth_map(16 * r + 1 : 24 * r, 1 : 1 * r) = 1; 

    Ground_truth_map(16*r + 1: 17*r, 1 : 13*r) = 1; 
    Ground_truth_map(3*r+1: 17*r, 12*r+1:13*r) = 1;
    Ground_truth_map(3*r+1: 4*r, 12*r+1: 33*r) = 1;
    Ground_truth_map(3*r+1: 24*r, 32*r+1: 33*r) = 1;
    Ground_truth_map(23*r+1: 24*r, 19*r+1 : 33*r) = 1;
    Ground_truth_map(23*r+1: 32*r, 19*r+1: 20*r) = 1;
    Ground_truth_map(31*r+1: 32*r, 12*r+1: 20*r) = 1;
    Ground_truth_map(23*r+1: 32*r, 12*r+1: 13*r) = 1;
    Ground_truth_map(23*r+1: 24*r, 1:13*r) = 1;

    Ground_truth_map(9*r+1: 10*r, 19*r+1: 24*r) = 1;
    Ground_truth_map(9*r+1: 15*r, 23*r+1: 24*r) = 1;
    Ground_truth_map(14*r+1: 15*r, 23*r+1: 28*r) = 1;
    Ground_truth_map(14*r+1: 17*r, 27*r+1: 28*r) = 1;
    Ground_truth_map(16*r+1: 17*r, 19*r+1: 28*r) = 1;
    Ground_truth_map(9*r+1: 17*r,19*r+1: 20*r) = 1;

    Ground_truth_map(9*r+1: 12*r, 25*r+1: 28*r) = 1;

    Ground_truth_map(18*r+1: 19*r, 2*r+1: 4*r) = 1;
    Ground_truth_map(18*r+1: 19*r, 6*r+1: 8*r) = 1;
    Ground_truth_map(18*r+1: 19*r,9*r+1: 12*r) = 1;

    Ground_truth_map(5*r+1: 8*r, 14*r+1: 15*r) = 1;
    Ground_truth_map(9*r+1: 11*r, 14*r+1: 15*r) = 1;

    Ground_truth_map(21*r+1: 22*r, 21*r+1: 24*r) = 1;
    Ground_truth_map(21*r+1: 22*r, 25*r+1: 27*r) = 1;
    Ground_truth_map(21*r+1: 22*r, 28*r+1: 31*r) = 1;

    Ground_truth_map(4*r+1: 5*r, 25*r+1: 28*r) = 1;
    Ground_truth_map(9*r+1: 10*r, 31*r+1: 32*r) = 1;
    Ground_truth_map(16*r+1: 17*r, 31*r+1: 32*r) = 1;
end
if FLAG == 1
    %% wall
    wallX = 8;
    Ground_truth_map(3*r+1: 7*r, (wallX-1)*r+1: wallX*r) = 1;
end
end