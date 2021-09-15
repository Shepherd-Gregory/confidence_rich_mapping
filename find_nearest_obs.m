function [ind_x, ind_y, uz_ray, uz_occ, uz_free, meas_range] = find_nearest_obs(map, ro_pose, ray_angle, param)
    z_ray = [];
    ind_x = [];
    ind_y = [];
    z_occ = [];
    rx = ro_pose(1);
    ry = ro_pose(2);
    rtheta = ro_pose(3);
    dz = 0.1 * 1/param.resol;
    z = 0;
    while z <= param.maxrange + dz
        x_z = z * cos(-ray_angle + rtheta) + rx ;
        y_z = z * sin(-ray_angle + rtheta) + ry ;
        ix_z = ceil(x_z * param.resol) + param.origin(1) ;
        iy_z = ceil(y_z * param.resol) + param.origin(2) ;
        ind_x = [ind_x; ix_z];
        ind_y = [ind_y; iy_z];
        % if hit the border
        if ix_z > param.mapSize * param.resol || ix_z < 1 || iy_z > param.mapSize * param.resol || iy_z < 1
            meas_range = z;
            break
        end
        % if not hit the border,check this cell
        ind_z = sub2ind(size(map), param.mapSize * param.resol + 1 - iy_z, ix_z); 
        z_ray = [z_ray; ind_z];
 
        if map(ind_z) == 1 % stop when hit 
%             meas_range = (1/param.mapResol) * norm([ind_ro(1), ind_ro(2)] - [ix_z, iy_z]);
            meas_range = z;
            z_occ = [z_occ; ind_z];
            break
        end
        % MAXRANGE
        if abs(z - param.maxrange) < dz
            meas_range = z;
            break
        end
       
        z = z + dz;
    end
    [uz_ray, mm] = unique(z_ray);
    [~, mm] = sort(mm);
    uz_ray = uz_ray(mm, :);
    [uz_occ, mmm] = unique(z_occ);
    [~, mmm] = sort(mmm);
    uz_occ = uz_occ(mmm, :);
    
    uz_free = uz_ray;
    if ~isempty(z_occ)
        for i = 1:length(z_occ)
            temp = find(uz_free == uz_occ(i));
            uz_free(temp) = [];
        end
        
    end
end