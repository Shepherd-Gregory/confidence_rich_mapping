function [myMap, com, bin_m, bel_m, exp_m ]  = confident_rich_mapping(V_ray, ray, myMap, ix_robot, iy_robot, com, bin_m, bel_m, exp_m, param)
    %% CRM
    res_m = param.res_m;
    vec_com = 1/res_m * (0:1:param.res_m-1);
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
        exp_m(V_ray(j)) = (1/param.res_m) * sum(vec_com.*bel_m(:,V_ray(j))');
    end
    if ray >= param.maxrange % maxrange case
        prod = 1;
        for j = 1:length(V_ray)
            prod = prod * (1- exp_m(V_ray(j)) / param.resol);
            f = prod;
            pz = 1/ param.maxrange;
            % compute SCM p(c_k| z{0:k},x{0:k})
            pc = pz * f;
            vec_scm(j) = pc;
            eta = eta + pc;
        end
        vec_scm = insertSort(vec_scm, min(vec_scm));
        eta = eta + 1/ param.maxrange;
    else % non-maxrange case
        for j = 1:length(V_ray)
            if j==1
                prod = 1;
            else
                prod = prod * (1- exp_m(V_ray(j-1)) / param.resol);
            end

            f = prod * exp_m(V_ray(j));
            [iy,ix] = ind2sub(param.size, V_ray(j));
            dist = (1/ param.resol)*norm([ix_robot, param.mapSize * param.resol + 1 - iy_robot] - [ix,iy]);
            pz = beamSensorModel(dist, ray, param);
            % compute SCM p(c_k| z{0:k},x{0:k})
            pc = pz * f;
            vec_scm(j) = pc;
            eta = eta + pc;
        end
    end

    % normalization
    vec_scm = vec_scm / eta;
    %% compute alpha and beta for belief
	beta = zeros(1,length(V_ray));
    alpha = zeros(1,length(V_ray));
    for j = 1:length(V_ray)
        beta(j) = 1 + exp_m(V_ray(j))/(1 - exp_m(V_ray(j)))*sum(vec_scm(j + 1:length(V_ray))) - vec_scm(j);
        alpha(j) = (1  - beta(j))/exp_m(V_ray(j));
    end

    % update the discretized belief
    for j = 1:size(V_ray)
        bel_m(:,V_ray(j)) = bel_m(:,V_ray(j)) .* (alpha(j) * vec_com'  +  beta(j) * ones(param.res_m,1));
		% normalization
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