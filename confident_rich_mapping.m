function [myMap, com, bin_m, bel_m, exp_m ]  = confident_rich_mapping(V_ray, ray, myMap, ix_robot, iy_robot, com, bin_m, bel_m, exp_m, param)
    %% CRM
    res_m = param.res_m;
    vec_com = 1/res_m * (1:1:param.res_m);
    eta = 0;
    vec_scm = zeros(length(V_ray));
    for j = 1:length(V_ray)
        com(V_ray(j)) = 2^myMap(V_ray(j))/(1 + 2^myMap(V_ray(j)));
        % update the binary map
        if com(V_ray(j)) >= 0.7
            bin_m(V_ray(j)) = 1;
        elseif com(V_ray(j)) <= 0.3
            bin_m(V_ray(j)) = -1;
        else
            bin_m(V_ray(j)) = 0;
        end
        % use a low resulution
        % 网格分辨率对概率的影响
        exp_m(V_ray(j)) = (1/param.res_m) * sum(vec_com.*bel_m(:,V_ray(j))');
    end
    if ray >= param.maxrange % maxrange 情况
        prod = 1;
        for j = 1:length(V_ray)
            % prod = prod * (1-exp_m(V_ray(j)));
            prod = prod * (1- exp_m(V_ray(j)) / param.resol);
            f = prod;
            % 传感模型计算每个栅格的p(z_k|c_k=m_i,x_k)
            % 考虑超出最大范围的概率
            pz = 1/ param.maxrange;
            % 计算SCM p(c_k| z{0:k},x{0:k})
            pc = pz * f;
            vec_scm(j) = pc;
            eta = eta + pc;
        end
        % 在maxrange情况下，添上界外的ck的概率
        vec_scm = insertSort(vec_scm, min(vec_scm));
        eta = eta + 1/ param.maxrange;
    else % 非 maxrange情况
        for j = 1:length(V_ray)
            if j==1
                prod = 1;
            else
                %                 prod = prod * (1-exp_m(V_ray(j-1)));
                prod = prod * (1- exp_m(V_ray(j-1)) / param.resol);
            end

            f = prod * exp_m(V_ray(j));

            % 计算地图上机器人和ray上各栅格之间的两点距离，
            % 由于网格对角线较长，该方法误差较大，造成鬼影的原因之一
            [iy,ix] = ind2sub(param.size, V_ray(j));
            dist = (1/ param.resol)*norm([ix_robot, param.mapSize * param.resol + 1 - iy_robot] - [ix,iy]);
            % 传感模型计算每个栅格的p(z_k|c_k=m_i,x_k)
            % 考虑超出最大范围的概率
            pz = beamSensorModel(dist, ray, param);
            % 计算SCM p(c_k| z{0:k},x{0:k})
            pc = pz * f;
            vec_scm(j) = pc;
            eta = eta + pc;
        end
    end

    % 归一化
    vec_scm = vec_scm / eta;
    %% 计算似然参数 alpha 和 beta
	beta = zeros(1,length(V_ray));
    alpha = zeros(1,length(V_ray));
    for j = 1:length(V_ray)
        beta(j) = 1 + exp_m(V_ray(j))/(1 - exp_m(V_ray(j)))*sum(vec_scm(j + 1:length(V_ray))) - vec_scm(j);
        % 对beta归一化？
		% beta(j) = max(beta(j), 1.01);
        alpha(j) = (1  - beta(j))/exp_m(V_ray(j));
    end

    % 计算更新信度值 ,此处乘上mi
    for j = 1:size(V_ray)
        bel_m(:,V_ray(j)) = bel_m(:,V_ray(j)) .* (alpha(j) * vec_com'  +  beta(j) * ones(param.res_m,1));
        % 归一化处理
        bel_m(:,V_ray(j)) = bel_m(:,V_ray(j)) * param.res_m / sum(bel_m(:,V_ray(j)));
    end
    %% DEBUG
    if size(find(ismissing(bel_m))) > 0
        pause;
    end

end

function y=insertSort(X, x)
    if x == X(1)
        y = [x; X];
    else
        y = [X; x];
    end
end