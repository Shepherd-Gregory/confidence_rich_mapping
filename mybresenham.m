function [ind_x, ind_y] = mybresenham(startPts, endPts)
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
coords = [];

dx = abs(endPts(1) - startPts(1));
dy = abs(endPts(2) - startPts(2));
x = startPts(1);
y = startPts(2);
if startPts(1) > endPts(1)
    sx = -1;
else
    sx = 1;
end
if startPts(2) > endPts(2)
    sy = -1;
else
    sy = 1;
end

if dx > dy
    err = dx / 2;
    while x ~= endPts(1)
        coords = [coords; x, y];
        err = err - dy;
        if err < 0
            y = y + sy;
            err = err + dx;
        end
        x = x + sx;
    end
else
    err = dy / 2;
    while y ~= endPts(2)
        coords = [coords; x, y];
        err = err - dx;
        if err < 0
            x = x + sx;
            err = err + dy;
        end
        y = y + sy;
    end
end
coords = [coords; x, y];
ind_x = coords(:,1);
ind_y = coords(:,2);
end