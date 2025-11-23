%% Base Performance

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_base_performance.mat";

load_system(model)

params = struct();
params.ref_angle = 2;

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 72; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 5;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";


out_table = zeros(eval_params.nbEval, eval_params.nbOut);
        
nb_done = 0;
while nb_done < eval_params.nbEval
    nb_remaining = eval_params.nbEval-nb_done;
    simins(1:nb_remaining) = Simulink.SimulationInput(model);
    for j = 1:length(simins)
        simins(j) = setVariable(simins(j), 'ref_angle', params.ref_angle, "Workspace", model);
        simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
        simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
        simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
    end
    
    out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
    
    for j = 1:nb_remaining
        if ~isnan(out(j).val(1)) % Simulation crashed -> relaunch
            nb_done = nb_done + 1;
            out_table(nb_done,:) = out(j).val;
        end
    end
end

save(save_path, "out_table", "params", "eval_params", "model");

%% Setup finer mismatch

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_montecarlo_finer.mat";

load_system(model)

params = struct();
params.ref_angle = 2;

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 72; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.5;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";


% Distributions
stddev_ratio = 2;
noise_names = {"amp_log_mismatch", "bias_mismatch", "gain_log_mismatch", "all_mismatch"};
noise_names_link = {["gain"], ["bias"], ["slope"], ["gain", "mod_gain", "bias", "timescale", "slope"]};

try
    old = load(save_path);
    
    noise_vec = old.noise_vec;
    start_k = old.k;
    start_i = old.i+1;
    out_table = old.out_table;
    out_data = old.out_data;
catch E
    disp("Old save incompatible, data will be erased");

    noise_vec = 0.0075:0.0075:0.075;
    start_k = 1;
    start_i = 1;
    out_table = zeros(length(noise_names), length(noise_vec), eval_params.nbEval, eval_params.nbOut);
    out_data = cell(length(noise_names), length(noise_vec), eval_params.nbEval);
end

%% Setup coarser mismatch

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_montecarlo_coarser.mat";

load_system(model)

params = struct();
params.ref_angle = 2;

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 72; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.5;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";


% Distributions
stddev_ratio = 2;
noise_names = {"amp_mod_log_mismatch", "timescale_log_mismatch"};
noise_names_link = {["mod_gain"], ["timescale"], ["slope"]};

try
    old = load(save_path);
    
    noise_vec = old.noise_vec;
    start_k = old.k;
    start_i = old.i+1;
    out_table = old.out_table;
    out_data = old.out_data;
catch E
    disp("Old save incompatible, data will be erased");

    noise_vec = 0.025:0.025:0.25;
    start_k = 1;
    start_i = 1;
    out_table = zeros(length(noise_names), length(noise_vec), eval_params.nbEval, eval_params.nbOut);
    out_data = cell(length(noise_names), length(noise_vec), eval_params.nbEval);
end

%% Mismatch Performance

import NeurCont.mismatch.apply_to_system
import NeurCont.mismatch.base_options

specific_mismatch_func = containers.Map( ...
    {'TransferFcn'},...
    { ...
        containers.Map(...
            {'timescale'},...
            { ...
                @(blPath, simins, model_work, mismatch_param) mismatch_dirty_diff_timescale(blPath, simins, model_work, mismatch_param.width, mismatch_param.std)...  
            }...
        )
    }...
);

for k = start_k:length(noise_names)
    name = noise_names_link{k};
    for i = start_i:length(noise_vec)
        noise_strength = noise_vec(i);
    
        mismatch_width = 2*noise_strength*stddev_ratio;
        mismatch_std = noise_strength;

        disp(noise_names{k} + ": noise = "+ noise_strength)
        
        nb_done = 0;
        while nb_done < eval_params.nbEval
            nb_remaining = eval_params.nbEval-nb_done;

            mismatch_param.width = 2*noise_strength*stddev_ratio;
            mismatch_param.std = noise_strength;

            simins = apply_to_system(model, mismatch_param, nb_remaining, base_options(), ...
                'mismatchFunctions', specific_mismatch_func,...
                'mismatchIncludeList',name, 'blockTypeIncludeList', 'all');
            for j = 1:length(simins)
                simins(j) = setVariable(simins(j), 'ref_angle', params.ref_angle, "Workspace", model);
                simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
                simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
                simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
            end
            
            out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
            
            for j = 1:nb_remaining
                if ~isnan(out(j).val(1)) % Simulation crashed -> relaunch
                    nb_done = nb_done + 1;
                    out_table(k,i,nb_done,:) = out(j).val;
                    out_data{k,i,nb_done} = out(j).data;
                end
            end
        end
    
        save(save_path, "noise_names", "noise_vec", "out_table", "out_data", "params", "eval_params", "model", "k", "i");
    end
    start_i = 1;
end

save(save_path, "noise_names", "noise_vec", "out_table", "out_data", "params", "eval_params", "model");

%% Add base perf to mismatch perf

% Choose source
perf_source = "full_mod_simplified_montecarlo_finer.mat";
%perf_source = "full_mod_simplified_montecarlo_coarser.mat";

base_perf = load("full_mod_simplified_base_performance.mat");
perf = load(perf_source);

out_table = zeros(size(perf.out_table, 1), size(perf.out_table, 2)+1, size(perf.out_table, 3), size(perf.out_table, 4));
out_table(:,2:size(perf.out_table, 2)+1,:,:) = perf.out_table;

for i = 1:size(perf.out_table, 1)
    out_table(i,1,:,:) = base_perf.out_table;
end

params = perf.params;
eval_params = perf.eval_params;
model = perf.model;
noise_vec = [0, perf.noise_vec];
noise_names = perf.noise_names;

save(perf_source, "noise_names", "noise_vec","out_table", "params", "eval_params", "model");

%% Mismatch Plots + Save to CSV

% Choose source
load("full_mod_simplified_montecarlo_finer.mat")
%load("full_mod_simplified_montecarlo_coarser.mat")

num_names = length(noise_names);

window_size = [500 1000];

for i = 1:num_names

    lims = {[-4, -1], [1, 4], [0, pi], [0.0, 0.3]};
    clamp = {[-15, 15], [0, 15], [], []};
    measure_names = {'g_{syn}', 'g_{s-}', 'Mean Absolute Amp Error', 'Mean Absolute Slip Error'};
    measure_file_names = {"gsyn", "gsm", "amp_error", "freq_error"};

    data_var = {squeeze(cat(3, out_table(i,:,:,1) , out_table(i,:,:,2))), squeeze(cat(3, out_table(i,:,:,5) , out_table(i,:,:,6))), squeeze(out_table(i,:,:,11)), squeeze(out_table(i,:,:,15))};

    f = figure;
    f.Position = [50 50 window_size];
    t = tiledlayout(length(data_var),1);
    
    axs = [];

    for j = 1:length(data_var)
        data = data_var{j};
        [average_std, average_mean] = std(data,0,2,"omitnan");
        quant = quantile(data,[0.05, 0.25, 0.5, 0.75, 0.95],2);
        bot_05_tab = quant(:,1);
        first_quartile_tab = quant(:,2);
        median_tab = quant(:,3);
        third_quartile_tab = quant(:,4);
        top_05_tab = quant(:,5);
        axs(j) = nexttile;
        hold on
        plot(axs(j), noise_vec, median_tab, "Color",ulgColors.BlueGreen.hex)
        scatter(axs(j), noise_vec, data, 5, ulgColors.BlueViolet.rgb, "filled")
        fill(axs(j), [noise_vec, flip(noise_vec)], [bot_05_tab', flip(top_05_tab')], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.1, "EdgeColor","none")
        fill(axs(j), [noise_vec, flip(noise_vec)], [first_quartile_tab', flip(third_quartile_tab')], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.3, "EdgeColor","none")
        %ylim(axs(j),lims{j})
        ylabel(axs(j),measure_names{j}, 'Interpreter', 'tex')
        

        if ~isempty(clamp{j})
            range = clamp{j};
            min_val = range(1);
            max_val = range(2);
            data(data < min_val) = NaN;
            data(data > max_val) = NaN;
        end

        T = array2table([noise_vec', average_mean, average_std, median_tab, bot_05_tab, first_quartile_tab, third_quartile_tab, top_05_tab], 'VariableNames', ["prob", "mean", "std", "median", "perc05", "perc25", "perc75", "perc95"]);
        writetable(T, "controller_mismatch_"+noise_names{i}+"_"+measure_file_names{j}+"_avg.csv");
        
        T = array2table([repelem(noise_vec,size(data,2))',reshape(data', 1,[])'], 'VariableNames', ["prob", "val"]);
        writetable(T, "controller_mismatch_"+noise_names{i}+"_"+measure_file_names{j}+"_pts.csv");
    end
    title(t,noise_names{i}, 'Interpreter', 'none')

    if true
        set(axs(3), 'YScale', 'log')
        set(axs(4), 'YScale', 'log')
    end
end

%% Event Loss Performances

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_loss.mat";

load_system(model)

params = struct();
params.ref_angle = 2;

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 72; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.5;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";


try
    old = load(save_path);
    
    prob_vec = old.prob_vec;
    start_i = old.i+1;
    out_table = old.out_table;
catch E
    disp("Old save incompatible, data erased");

    prob_vec = 0:0.05:0.75;
    start_i = 1;
    out_table = zeros(length(prob_vec), eval_params.nbEval, eval_params.nbOut);
end


for i = start_i:length(prob_vec)
    loss_prob = prob_vec(i);
    
    nb_done = 0;
    while nb_done <= eval_params.nbEval
        simins(1:eval_params.nbEval-nb_done) = Simulink.SimulationInput(model);
        for j = 1:length(simins)
            simins(j) = setVariable(simins(j), 'prob_event_loss', loss_prob, "Workspace", model);
            simins(j) = setVariable(simins(j), 'ref_angle', params.ref_angle, "Workspace", model);
            simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
            simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
            simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
        end
        
        out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
        
        for j = 1:eval_params.nbEval
            if ~isnan(out(j).val(1))
                nb_done = nb_done + 1;
                out_table(i,nb_done,:) = out(j).val;
            end
        end
    end

    save(save_path, "prob_vec","out_table", "params", "eval_params", "model", "i");
end

save(save_path, "prob_vec","out_table", "params", "eval_params", "model");

%% Event Loss Plots + Save to CSV

save_path = "full_mod_simplified_loss.mat";
load(save_path)

lims = {[0, pi], [0.0, 1], [-3, -1], [2, 5]};
measure_names = {'Mean Absolute Amp Error', 'Mean Absolute Slip Error', 'g_{syn}', 'g_{s-}'};
measure_file_names = {"amp_error", "freq_error", "gsyn", "gsm"};

data_var = {squeeze(out_table(:,:,11)), squeeze(out_table(:,:,15)), squeeze(cat(2, out_table(:,:,1) , out_table(:,:,2))), squeeze(cat(2, out_table(:,:,5) , out_table(:,:,6)))};


window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];
t = tiledlayout(length(data_var),1);
for j = 1:length(data_var)
    axs(j) = nexttile;
    hold on

    data = data_var{j};
    [average_std, average_mean] = std(data,0,2,"omitnan");
    quant = quantile(data,[0.05, 0.25, 0.5, 0.75, 0.95],2)';
    bot_05_tab = quant(1,:);
    first_quartile_tab = quant(2,:);
    median_tab = quant(3,:);
    third_quartile_tab = quant(4,:);
    top_05_tab = quant(5,:);
    hold on
    plot(prob_vec, median_tab, "Color",ulgColors.BlueGreen.hex)
    scatter(prob_vec, data', 5, ulgColors.BlueViolet.rgb, "filled")
    fill([prob_vec, flip(prob_vec)], [bot_05_tab, flip(top_05_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.1, "EdgeColor","none")
    fill([prob_vec, flip(prob_vec)], [first_quartile_tab, flip(third_quartile_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.3, "EdgeColor","none")
    ylim(lims{j})
    ylabel(measure_names{j}, 'Interpreter', 'tex')

    T = array2table([prob_vec', average_mean, average_std, median_tab', bot_05_tab', first_quartile_tab', third_quartile_tab', top_05_tab'], 'VariableNames', ["prob", "mean", "std", "median", "perc05", "perc25", "perc75", "perc95"]);
    writetable(T, "event_loss_"+measure_file_names{j}+"_avg.csv");
    
    T = array2table([repelem(prob_vec,size(data,2))',reshape(data, 1,[])'], 'VariableNames', ["prob", "val"]);
    writetable(T, "event_loss_"+measure_file_names{j}+"_pts.csv");
end

if true
    set(axs(1), 'YScale', 'log')
    set(axs(2), 'YScale', 'log')
end


%% Angle Reference Performance

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_ref_angle.mat";

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 24; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.5;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";
eval_params.save_path = "param_chart_save_tmp.mat";

params = struct();

var_params = struct();
var_params.ref_angle = 0.1:0.1:3.1;

out_table = parameterChart(params, var_params, model, eval_params, 72);

save(save_path, "out_table", "params", "var_params", "eval_params", "model");

%% Angle Reference Plots + Save to CSV

save_path = "full_mod_simplified_ref_angle.mat";
load(save_path)


[~, max_num_rows_nan] = cleanUpResultArray(out_table, 3, 2, 1, 1, "onlyNaN");
[~, max_num_rows_out] = cleanUpResultArray(out_table, 3, 2, 1, 1, "outliers");

lims = {[0, pi], [0.0, 1], [-4, -1], [1.5, 5]};
measure_names = {'Mean Absolute Amp Error', 'Mean Absolute Slip Error', 'g_{syn}', 'g_{s-}'};
measure_file_names = {"amp_error", "freq_error", "gsyn", "gsm"};

data_var = {squeeze(out_table(11,:,:)), squeeze(out_table(15,:,:)), squeeze(cat(2, out_table(1,:,:) , out_table(2,:,:))), squeeze(cat(2, out_table(5,:,:) , out_table(6,:,:)))};


window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];
t = tiledlayout(length(data_var),1);
for j = 1:length(data_var)
    axs(j) = nexttile;
    hold on

    data = data_var{j};
    [average_std, average_mean] = std(data,0,1,"omitnan");
    quant = quantile(data,[0.05, 0.25, 0.5, 0.75, 0.95],1);
    bot_05_tab = quant(1,:);
    first_quartile_tab = quant(2,:);
    median_tab = quant(3,:);
    third_quartile_tab = quant(4,:);
    top_05_tab = quant(5,:);
    hold on
    plot(var_params.ref_angle, median_tab, "Color",ulgColors.BlueGreen.hex)
    scatter(var_params.ref_angle, data', 5, ulgColors.BlueViolet.rgb, "filled")
    fill([var_params.ref_angle, flip(var_params.ref_angle)], [bot_05_tab, flip(top_05_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.1, "EdgeColor","none")
    fill([var_params.ref_angle, flip(var_params.ref_angle)], [first_quartile_tab, flip(third_quartile_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.3, "EdgeColor","none")
    ylim(lims{j})
    ylabel(measure_names{j}, 'Interpreter', 'tex')

    T = array2table([var_params.ref_angle', average_mean', average_std', median_tab', bot_05_tab', first_quartile_tab', third_quartile_tab', top_05_tab'], 'VariableNames', ["prob", "mean", "std", "median", "perc05", "perc25", "perc75", "perc95"]);
    writetable(T, "ref_angle_"+measure_file_names{j}+"_avg.csv");
    
    T = array2table([repelem(var_params.ref_angle,size(data,1))',reshape(data, 1,[])'], 'VariableNames', ["prob", "val"]);
    writetable(T, "ref_angle_"+measure_file_names{j}+"_pts.csv");
end

if true
    set(axs(1), 'YScale', 'log')
    set(axs(2), 'YScale', 'log')
end



%% Mod Strength Performance

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_mod_gains.mat";

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 24; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.55;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";
eval_params.save_path = "param_chart_save_tmp.mat";

params = struct();
params.ref_angle = 2;

var_params = struct();
var_params.mod_gain = 0.2:0.2:3.8;

out_table = parameterChartLinked(params, var_params, model, @link_mod_gains, eval_params, 72);

save(save_path, "out_table", "params", "var_params", "eval_params", "model");

%% Mod Strength Plots + Save to CSV

save_path = "full_mod_simplified_mod_gains.mat";
load(save_path)

lims = {[0, 0.1], [0.0, 0.1], [-3, -2], [2, 3]};
measure_names = {'Mean Absolute Amp Error', 'Mean Absolute Slip Error', 'g_{syn}', 'g_{s-}'};
measure_file_names = {"amp_error", "freq_error", "gsyn", "gsm"};

data_var = {squeeze(out_table(11,:,:)), squeeze(out_table(15,:,:)), squeeze(cat(2, out_table(1,:,:) , out_table(2,:,:))), squeeze(cat(2, out_table(5,:,:) , out_table(6,:,:)))};


window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];
t = tiledlayout(length(data_var),1);
for j = 1:length(data_var)
    axs(j) = nexttile;
    hold on

    data = data_var{j};
    [average_std, average_mean] = std(data,0,1,"omitnan");
    quant = quantile(data,[0.05, 0.25, 0.5, 0.75, 0.95],1);
    bot_05_tab = quant(1,:);
    first_quartile_tab = quant(2,:);
    median_tab = quant(3,:);
    third_quartile_tab = quant(4,:);
    top_05_tab = quant(5,:);
    hold on
    plot(var_params.mod_gain, median_tab, "Color",ulgColors.BlueGreen.hex)
    scatter(var_params.mod_gain, data', 5, ulgColors.BlueViolet.rgb, "filled")
    fill([var_params.mod_gain, flip(var_params.mod_gain)], [bot_05_tab, flip(top_05_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.1, "EdgeColor","none")
    fill([var_params.mod_gain, flip(var_params.mod_gain)], [first_quartile_tab, flip(third_quartile_tab)], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.3, "EdgeColor","none")
    ylim(lims{j})
    ylabel(measure_names{j}, 'Interpreter', 'tex')

    T = array2table([var_params.mod_gain', average_mean', average_std', median_tab', bot_05_tab', first_quartile_tab', third_quartile_tab', top_05_tab'], 'VariableNames', ["prob", "mean", "std", "median", "perc05", "perc25", "perc75", "perc95"]);
    writetable(T, "mod_gains_"+measure_file_names{j}+"_avg.csv");
    
    T = array2table([repelem(var_params.mod_gain,size(data,1))',reshape(data, 1,[])'], 'VariableNames', ["prob", "val"]);
    writetable(T, "mod_gains_"+measure_file_names{j}+"_pts.csv");
end

if true
    set(axs(1), 'YScale', 'log')
    set(axs(2), 'YScale', 'log')
end


%% Plant Params Performance

model = "F_mod_A_mod_gsyn_top_trig_mism_sub_simplified";
save_path = "full_mod_simplified_plant_params.mat";

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 24; % Number of Cores Available = 12
eval_params.StopTime = 600;
eval_params.func = @postsim_full_mod_metrics;
eval_params.nbOut = 16;
% Computation params
eval_params.StartAnalyseTime = 400;
eval_params.peakProminence = 0.1;
eval_params.eventTresh = 0.55;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";
eval_params.save_path = "param_chart_save_tmp.mat";

params = struct();
params.ref_angle = 2;


value_names = {'mass_density', 'pendulum_length', 'friction_factor'};
value_base = {951, 0.4267, 1};

eval_values = logspace(log10(0.5), log10(2),15);

out_table = zeros(length(value_names), eval_params.nbOut, eval_params.nbEval, length(eval_values));

for i = 1:length(value_names)
    name = value_names{i};
    value = value_base{i};

    var_params = struct();
    var_params.(name) = value*eval_values;
    
    out_table(i,:,:,:) = parameterChart(params, var_params, model, eval_params, 72);
end

save(save_path, "out_table", "params", "eval_params", "model", "value_names", "value_base", "eval_values");


%% Plant Params Plots + Save to CSV

% vertical plots
load("full_mod_simplified_plant_params.mat")

num_names = length(value_names);

window_size = [500 1000];

for i = 1:num_names

    lims = {[-4, -1], [1, 4], [0, pi], [0.0, 0.3]};
    clamp = {[-15, 15], [0, 15], [], []};
    measure_names = {'g_{syn}', 'g_{s-}', 'Mean Absolute Amp Error', 'Mean Absolute Slip Error'};
    measure_file_names = {"gsyn", "gsm", "amp_error", "freq_error"};

    data_var = {squeeze(cat(3, out_table(i,1,:,:) , out_table(i,2,:,:))), squeeze(cat(3, out_table(i,5,:,:) , out_table(i,5,:,:))), squeeze(out_table(i,11,:,:)), squeeze(out_table(i,15,:,:))};
    
    switch value_names{i}
        case 'mass_density'
            radius = 0.05;
            value_vec = value_base{i}*eval_values*4*pi*(radius^3)/3;
        case 'friction_factor'
            deg_sec_factor = 5e-04;
            value_vec = value_base{i}*eval_values*deg_sec_factor*180/pi;
        otherwise
            value_vec = value_base{i}*eval_values;
    end

    f = figure;
    f.Position = [50 50 window_size];
    t = tiledlayout(length(data_var),1);
    
    axs = [];

    for j = 1:length(data_var)
        data = data_var{j}';
        [average_std, average_mean] = std(data,0,2,"omitnan");
        quant = quantile(data,[0.05, 0.25, 0.5, 0.75, 0.95],2);
        bot_05_tab = quant(:,1);
        first_quartile_tab = quant(:,2);
        median_tab = quant(:,3);
        third_quartile_tab = quant(:,4);
        top_05_tab = quant(:,5);
        axs(j) = nexttile;
        hold on
        plot(axs(j), value_vec, median_tab, "Color",ulgColors.BlueGreen.hex)
        scatter(axs(j), value_vec, data, 5, ulgColors.BlueViolet.rgb, "filled")
        fill(axs(j), [value_vec, flip(value_vec)], [bot_05_tab', flip(top_05_tab')], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.1, "EdgeColor","none")
        fill(axs(j), [value_vec, flip(value_vec)], [first_quartile_tab', flip(third_quartile_tab')], ulgColors.BlueGreenLight.rgb, "FaceAlpha", 0.3, "EdgeColor","none")
        ylim(axs(j),lims{j})
        ylabel(axs(j),measure_names{j}, 'Interpreter', 'tex')
        

        if ~isempty(clamp{j})
            range = clamp{j};
            min_val = range(1);
            max_val = range(2);
            data(data < min_val) = NaN;
            data(data > max_val) = NaN;
        end

        T = array2table([value_vec', average_mean, average_std, median_tab, bot_05_tab, first_quartile_tab, third_quartile_tab, top_05_tab], 'VariableNames', ["prob", "mean", "std", "median", "perc05", "perc25", "perc75", "perc95"]);
        writetable(T, "plant_variability_"+value_names{i}+"_"+measure_file_names{j}+"_avg.csv");
        
        T = array2table([repelem(value_vec,size(data,2))',reshape(data', 1,[])'], 'VariableNames', ["prob", "val"]);
        writetable(T, "plant_variability_"+value_names{i}+"_"+measure_file_names{j}+"_pts.csv");
    end
    title(t,value_names{i}, 'Interpreter', 'none')

    if true
        set(axs(3), 'YScale', 'log')
        set(axs(4), 'YScale', 'log')
    end
end



%% Script Local Functions

function simins = mismatch_dirty_diff_timescale(blPath, simins, model_work, mismatch_width, mismatch_std)
    name = 'Denominator';
    try
        value = model_work.evalin(get_param(blPath, name));
    catch E
        error(['Param value "' name '" from block "' blPath '" cannot be evaluated.']);
    end
    for j = 1:length(simins)
        new_value = monteLogCutNormal(value(1),mismatch_std,mismatch_width);
        simins(j) = simins(j).setBlockParameter(blPath, name, ['[' num2str(new_value, '%.5e') ', 1]']);
    end
end

function value = monteLogCutNormal(center, log_stddev, log_span)
    if 0 > log_stddev
        error("monteLogCUtNormal expected log_stddev >= 0");
    end
    if 0 > log_span
        error("monteLogCUtNormal expected log_span >= 0");
    end
    value = center * 10^distrib(log_stddev,log_span/2);
end

function r = distrib(stddev, cutoff)
	r = stddev*randn(); % Normally distributed ~ N(0,std^2)
    if r < -cutoff || r > cutoff
        r = (rand() - 0.5)*2*cutoff; % Uniformly distributed ~ U(-cutoff, cutoff)
    end
end

function cur_params = link_mod_gains(cur_params, tmp_params)
    gain = tmp_params.mod_gain;

    cur_params.mod_strength_gsm = gain;
    cur_params.mod_strength_gsyn = gain;
end