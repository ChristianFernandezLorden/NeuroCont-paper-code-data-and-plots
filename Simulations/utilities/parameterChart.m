function out_val = parameterChart(params, var_params, model, eval_params, nb_par)
    if strcmp(params, "resume")
        data = load(eval_params.save_path);
        out_val = data.out_val;
        nb_simed = data.nb_simed;
        total_sim = data.total_sim;
        pos_vec = data.pos_vec;
        ind = data.ind;
        pos = data.pos;
        dim = data.dim;
        params = data.params;
        var_params = data.var_params;
        model = data.model;
        eval_params = data.eval_params;

        load_system(model)
    
        fixed_param_fields = fieldnames(params);
        updatable_params_fields = fieldnames(var_params);
    else
        load_system(model)
    
        fixed_param_fields = fieldnames(params);
        updatable_params_fields = fieldnames(var_params);
        
        dim = [];
        for i = 1:length(updatable_params_fields)
            name = updatable_params_fields{i};
            dim(i) = length(var_params.(name));
        end
        
        out_dim = [eval_params.nbOut, eval_params.nbEval, dim];
        out_val = zeros(out_dim);
        
        ind = 1;
        nb_simed = 0;
        pos = ones(length(dim), 1);
        pos_vec = zeros(nb_par, 1);
    
        nb_points = prod(dim, 'all');
        total_sim = eval_params.nbEval*nb_points;
    end
    disp("Number of Sims Done : "+string(nb_simed)+"/"+string(total_sim))
    while nb_simed < total_sim
        nb_remaining = total_sim - nb_simed;
        nb_sim = min(nb_remaining, nb_par);

        sim_in(1:nb_sim) = Simulink.SimulationInput(model);
        sim_in = sim_in(1:nb_sim);
        for j = 1:nb_sim
            cur_params = struct();

            for k = 1:length(fixed_param_fields)
                name = fixed_param_fields{k};
                cur_params.(name) = params.(name);
            end
    
            for k = 1:length(updatable_params_fields)
                name = updatable_params_fields{k};
                vec_param = var_params.(name);
                cur_params.(name) = vec_param(pos(k));
            end

            sim_in(j) = load_params(sim_in(j), cur_params, model);
            
            if isfield(eval_params, 'simulation_mode')
                sim_in(j) = setModelParameter(sim_in(j), "SimulationMode", eval_params.simulation_mode);
            end
            sim_in(j) = setModelParameter(sim_in(j), "StopTime", num2str(eval_params.StopTime));
    

            sim_in(j) = setPostSimFcn(sim_in(j), @(o) (eval_params.func(o, eval_params, cur_params)));
            
            pos_vec(j) = ind;
            ind = ind+1;

            %{
            ind = ind+1;
            sim_in(ind) = sim_tmp;
    
            if (ind == nb_par)
                ind = 0;
                out = parsim(sim_in);
                disp(string(i)+"/"+string(nb_points))
                for k = 1:nb_par
                    out_val(:, pos_vec(k)) = reshape(out_val(:, pos_vec(k)), 1, []) + reshape(out(k).val, 1, []);
                end
            end
            %}

            if mod(ind,eval_params.nbEval) == 0
                p = 1;
                while p <= length(pos)
                    pos(p) = pos(p)+1;
                    if pos(p) > dim(p)
                        pos(p) = 1;
                        p = p+1;
                    else
                        break
                    end
                end
            end
        end

        out = simulate(sim_in, eval_params);
        for k = 1:nb_sim
            out_val(:, pos_vec(k)) = reshape(out_val(:, pos_vec(k)), 1, []) + reshape(out(k).val, 1, []);
        end
        nb_simed = nb_simed + nb_sim;
        disp("Number of Sims Done : "+string(nb_simed)+"/"+string(total_sim))

        if isfield(eval_params, 'save_path')
            save(eval_params.save_path, "out_val", "nb_simed", "total_sim",...
                 "pos_vec", "ind", "pos", "dim", ...
                 "params", "var_params", "model", "eval_params", "nb_par")
        end
    end

    %out_val = out_val/eval_params.nbEval;
end

function out = simulate(sim_in, eval_params)
    if isfield(eval_params, 'fast_restart')
        out = parsim(sim_in, "UseFastRestart", eval_params.fast_restart);
    else
        out = parsim(sim_in);
    end
end

function simin = load_params(simin, params, model)
    param_fields = fieldnames(params);
    for i = 1:length(param_fields)
        name = param_fields{i};
        simin = setVariable(simin, name, params.(name), "Workspace", model);
    end
end