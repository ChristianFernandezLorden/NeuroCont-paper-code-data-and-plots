%% Baseline Simulation (with a change of 'a' parameter)

model = "simplified_STG";

load_system(model)

t_span = [10, 15];

simin = Simulink.SimulationInput(model);
simin = simin.setVariable('a',0.01, "Workspace", model);
out = sim(simin);

values = squeeze(out.events.signals.values);
[red_one, red_time_one] = reduce_points(values(:,1), out.events.time, 0.001);
[red_two, red_time_two] = reduce_points(values(:,2), out.events.time, 0.001);
[red_three, red_time_three] = reduce_points(values(:,3), out.events.time, 0.001);

i_start = find(red_time_one> t_span(1),1)-1;
i_end = find(red_time_one> t_span(2),1);
red_time_one = red_time_one(i_start:i_end);
red_one = red_one(i_start:i_end);
red_spike_time_one = spike_times_finder(red_one, red_time_one);
red_spike_one = ones(size(red_spike_time_one));

i_start = find(red_time_two> t_span(1),1)-1;
i_end = find(red_time_two> t_span(2),1);
red_time_two = red_time_two(i_start:i_end);
red_two = red_two(i_start:i_end);
red_spike_time_two = spike_times_finder(red_two, red_time_two);
red_spike_two = ones(size(red_spike_time_two));

i_start = find(red_time_three> t_span(1),1)-1;
i_end = find(red_time_three> t_span(2),1);
red_time_three = red_time_three(i_start:i_end);
red_three = red_three(i_start:i_end);
red_spike_time_three = spike_times_finder(red_three, red_time_three);
red_spike_three = ones(size(red_spike_time_three));

figure;
hold on
stem(red_spike_time_one, red_spike_one, "Color",ulgColors.BlueGreen.hex);
stem(red_spike_time_two, red_spike_two, "Color",ulgColors.YellowOrange.hex);
stem(red_spike_time_three, red_spike_three, "Color",ulgColors.Red.hex);
xlim(t_span)

red_spike_time_one = red_spike_time_one - t_span(1);
red_spike_time_two = red_spike_time_two - t_span(1);
red_spike_time_three = red_spike_time_three - t_span(1);

T = array2table([red_spike_time_one, red_spike_one], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron1_a015.csv");
T = array2table([red_spike_time_two, red_spike_two], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron2_a015.csv");
T = array2table([red_spike_time_three, red_spike_three], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron3_a015.csv");


simin = Simulink.SimulationInput(model);
simin = simin.setVariable('a',0.3, "Workspace", model);
out = sim(simin);

values = squeeze(out.events.signals.values);
[red_one, red_time_one] = reduce_points(values(:,1), out.events.time, 0.001);
[red_two, red_time_two] = reduce_points(values(:,2), out.events.time, 0.001);
[red_three, red_time_three] = reduce_points(values(:,3), out.events.time, 0.001);

i_start = find(red_time_one> t_span(1),1)-1;
i_end = find(red_time_one> t_span(2),1);
red_time_one = red_time_one(i_start:i_end);
red_one = red_one(i_start:i_end);
red_spike_time_one = spike_times_finder(red_one, red_time_one);
red_spike_one = ones(size(red_spike_time_one));

i_start = find(red_time_two> t_span(1),1)-1;
i_end = find(red_time_two> t_span(2),1);
red_time_two = red_time_two(i_start:i_end);
red_two = red_two(i_start:i_end);
red_spike_time_two = spike_times_finder(red_two, red_time_two);
red_spike_two = ones(size(red_spike_time_two));

i_start = find(red_time_three> t_span(1),1)-1;
i_end = find(red_time_three> t_span(2),1);
red_time_three = red_time_three(i_start:i_end);
red_three = red_three(i_start:i_end);
red_spike_time_three = spike_times_finder(red_three, red_time_three);
red_spike_three = ones(size(red_spike_time_three));

figure;
hold on
stem(red_spike_time_one, red_spike_one, "Color",ulgColors.BlueGreen.hex);
stem(red_spike_time_two, red_spike_two, "Color",ulgColors.YellowOrange.hex);
stem(red_spike_time_three, red_spike_three, "Color",ulgColors.Red.hex);
xlim(t_span)

red_spike_time_one = red_spike_time_one - t_span(1);
red_spike_time_two = red_spike_time_two - t_span(1);
red_spike_time_three = red_spike_time_three - t_span(1);

T = array2table([red_spike_time_one, red_spike_one], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron1_a05.csv");
T = array2table([red_spike_time_two, red_spike_two], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron2_a05.csv");
T = array2table([red_spike_time_three, red_spike_three], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron3_a05.csv");

%% Baseline Performance

model = "simplified_STG";
save_path = "stg_baseline_neuron1_out.mat";

load_system(model)

params = struct();

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 144; % Number of Cores Available = 12
eval_params.StopTime = 60;
eval_params.func = @postsim_neuronal_signal_analysis;
eval_params.nbOut = 9;
% Computation params
eval_params.StartAnalyseTime = 20;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";

out_table = zeros(eval_params.nbEval, eval_params.nbOut);

out_height_vec = [];
out_burst_freq_vec = [];
out_spikeburst_vec = [];

simins(1:eval_params.nbEval) = Simulink.SimulationInput(model);

%simins = mismatch_system(model, mismatch_width, mismatch_std, eval_params.nbEval);
for j = 1:length(simins)
    simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
    simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
    simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
    simins(j) = setVariable(simins(j), 'a',0.01, "Workspace", model);
end

out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
for i = 1:eval_params.nbEval
    out_table(i,:) = out(i).val;
    out_height_vec = [out_height_vec; out(i).data.spike_heights];
    out_burst_freq_vec = [out_burst_freq_vec; 1./out(i).data.bursts_cycle];
    out_spikeburst_vec = [out_spikeburst_vec; out(i).data.bursts_nb_spikes];
end

save(save_path,"out_table", "out_height_vec", "out_burst_freq_vec", "out_spikeburst_vec", "params", "eval_params", "model");

%% Baseline Plots + Save to CSV

save_path = "stg_baseline_neuron1_out.mat";

load(save_path)

window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];
t = tiledlayout(4,1);

axs(1) = nexttile;
histogram(out_spikeburst_vec)
ylabel("nb spikes", 'Interpreter', 'tex')

axs(2) = nexttile;
histogram(out_burst_freq_vec)
ylabel("burst freq", 'Interpreter', 'tex')

axs(3) = nexttile;
stem(out_table(:,7), ones(eval_params.nbEval,1))
ylabel("burst length", 'Interpreter', 'tex')

axs(4) = nexttile;
histogram(out_height_vec)
ylabel("spike height", 'Interpreter', 'tex')


T = array2table([out_table(:,3), out_table(:,6), out_table(:,7), out_table(:,9), ones(eval_params.nbEval,1)], 'VariableNames', ["nbspikes", "freq", "length", "height", "val"]);
writetable(T, "stg_baseline_avg.csv");

T = array2table([out_spikeburst_vec, ones(length(out_spikeburst_vec),1)], 'VariableNames', ["nbspikes", "val"]);
writetable(T, "stg_baseline_nbspike.csv");
T = array2table([out_burst_freq_vec, ones(length(out_burst_freq_vec),1)], 'VariableNames', ["freq", "val"]);
writetable(T, "stg_baseline_burstfreq.csv");
T = array2table([out_height_vec, ones(length(out_height_vec),1)], 'VariableNames', ["height", "val"]);
writetable(T, "stg_baseline_spikeheight.csv");
 
%% Mismatch Performance

model = "simplified_STG";
save_path = "stg_mismatch_neuron1_out.mat";

load_system(model)

params = struct();

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 12;%144; % Number of Cores Available = 12
eval_params.StopTime = 60;
eval_params.func = @postsim_neuronal_signal_analysis;
eval_params.nbOut = 9;
% Computation params
eval_params.StartAnalyseTime = 20;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";

out_table = zeros(eval_params.nbEval, eval_params.nbOut);
out_height_vec = [];
out_burst_freq_vec = [];
out_spikeburst_vec = [];

stddev_ratio = 2;
noise_strength = 0.025;

mismatch_param.width = 2*noise_strength*stddev_ratio;
mismatch_param.std = noise_strength;

import NeuroCont.mismatch.apply_to_system
import NeuroCont.mismatch.base_options

simins = apply_to_system(model, mismatch_param, eval_params.nbEval, base_options(), 'mismatchIncludeList', 'all', 'blockTypeIncludeList', 'all');

for j = 1:length(simins)
    simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
    simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
    simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
    simins(j) = setVariable(simins(j), 'a',0.01, "Workspace", model);
end

out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
for i = 1:eval_params.nbEval
    out_table(i,:) = out(i).val;
    out_height_vec = [out_height_vec; out(i).data.spike_heights];
    out_burst_freq_vec = [out_burst_freq_vec; 1./out(i).data.bursts_cycle];
    out_spikeburst_vec = [out_spikeburst_vec; out(i).data.bursts_nb_spikes];
end

save(save_path, "noise_strength", "stddev_ratio","out_table", "out_height_vec", "out_burst_freq_vec", "out_spikeburst_vec", "params", "eval_params", "model");

%% Mismatch Plots + Save to CSV

save_path = "stg_mismatch_neuron1_out.mat";

load(save_path)

window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];
t = tiledlayout(4,1);

axs(1) = nexttile;
histogram(out_spikeburst_vec)
ylabel("nb spikes", 'Interpreter', 'tex')

axs(2) = nexttile;
histogram(out_burst_freq_vec)
ylabel("burst freq", 'Interpreter', 'tex')

axs(3) = nexttile;
stem(out_table(:,7), ones(eval_params.nbEval,1))
ylabel("burst length", 'Interpreter', 'tex')

axs(4) = nexttile;
histogram(out_height_vec)
ylabel("spike height", 'Interpreter', 'tex')


T = array2table([out_table(:,3), out_table(:,6), out_table(:,7), out_table(:,9), ones(eval_params.nbEval,1)], 'VariableNames', ["nbspikes", "freq", "length", "height", "val"]);
writetable(T, "stg_mismatch_avg.csv");

T = array2table([out_spikeburst_vec, ones(length(out_spikeburst_vec),1)], 'VariableNames', ["nbspikes", "val"]);
writetable(T, "stg_mismatch_nbspike.csv");
T = array2table([out_burst_freq_vec, ones(length(out_burst_freq_vec),1)], 'VariableNames', ["freq", "val"]);
writetable(T, "stg_mismatch_burstfreq.csv");
T = array2table([out_height_vec, ones(length(out_height_vec),1)], 'VariableNames', ["height", "val"]);
writetable(T, "stg_mismatch_spikeheight.csv");

%% Mismatch Single Simulation to CSV

model = "simplified_STG";

load_system(model)

t_span = [10, 15];

stddev_ratio = 2;
noise_strength = 0.025;

mismatch_param.width = 2*noise_strength*stddev_ratio;
mismatch_param.std = noise_strength;

import NeuroCont.mismatch.apply_to_system
import NeuroCont.mismatch.base_options

simin = apply_to_system(model, mismatch_param, 1, base_options(), 'mismatchIncludeList', 'all', 'blockTypeIncludeList', 'all');
simin = simin.setVariable('a',0.01, "Workspace", model);
out = sim(simin);

values = squeeze(out.events.signals.values);
[red_one, red_time_one] = reduce_points(values(:,1), out.events.time, 0.001);
[red_two, red_time_two] = reduce_points(values(:,2), out.events.time, 0.001);
[red_three, red_time_three] = reduce_points(values(:,3), out.events.time, 0.001);

i_start = find(red_time_one> t_span(1),1)-1;
i_end = find(red_time_one> t_span(2),1);
red_time_one = red_time_one(i_start:i_end);
red_one = red_one(i_start:i_end);
red_spike_time_one = spike_times_finder(red_one, red_time_one);
red_spike_one = ones(size(red_spike_time_one));

i_start = find(red_time_two> t_span(1),1)-1;
i_end = find(red_time_two> t_span(2),1);
red_time_two = red_time_two(i_start:i_end);
red_two = red_two(i_start:i_end);
red_spike_time_two = spike_times_finder(red_two, red_time_two);
red_spike_two = ones(size(red_spike_time_two));

i_start = find(red_time_three> t_span(1),1)-1;
i_end = find(red_time_three> t_span(2),1);
red_time_three = red_time_three(i_start:i_end);
red_three = red_three(i_start:i_end);
red_spike_time_three = spike_times_finder(red_three, red_time_three);
red_spike_three = ones(size(red_spike_time_three));

figure;
hold on
stem(red_spike_time_one, red_spike_one, "Color",ulgColors.BlueGreen.hex);
stem(red_spike_time_two, 2*red_spike_two, "Color",ulgColors.YellowOrange.hex);
stem(red_spike_time_three, 3*red_spike_three, "Color",ulgColors.Red.hex);
xlim(t_span)

red_spike_time_one = red_spike_time_one - t_span(1);
red_spike_time_two = red_spike_time_two - t_span(1);
red_spike_time_three = red_spike_time_three - t_span(1);

T = array2table([red_spike_time_one, red_spike_one], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron1_mismatch.csv");
T = array2table([red_spike_time_two, red_spike_two], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron2_mismatch.csv");
T = array2table([red_spike_time_three, red_spike_three], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron3_mismatch.csv");

%% Event Loss Performance

model = "simplified_STG";
save_path = "stg_loss_neuron1_out.mat";

load_system(model)

params = struct();

eval_params = struct();

% Montecarlo Params
eval_params.nbEval = 144; % Number of Cores Available = 12
eval_params.StopTime = 60;
eval_params.func = @postsim_neuronal_signal_analysis;
eval_params.nbOut = 9;
% Computation params
eval_params.StartAnalyseTime = 20;
eval_params.simulation_mode = "rapid-accelerator";
eval_params.fast_restart = "off";

out_table = zeros(eval_params.nbEval, eval_params.nbOut);
out_height_vec = [];
out_burst_freq_vec = [];
out_spikeburst_vec = [];

loss_prob = 0.2;

simins(1:eval_params.nbEval) = Simulink.SimulationInput(model);
for j = 1:length(simins)
    simins(j) = setVariable(simins(j), 'loss_prob', loss_prob, "Workspace", model);
    simins(j) = setVariable(simins(j), 'a',0.01, "Workspace", model);
    simins(j) = setModelParameter(simins(j), "SimulationMode", eval_params.simulation_mode);
    simins(j) = setModelParameter(simins(j), "StopTime", num2str(eval_params.StopTime));
    simins(j) = setPostSimFcn(simins(j), @(o) (eval_params.func(o, eval_params, params)));
end

out = parsim(simins, "UseFastRestart", eval_params.fast_restart);
for i = 1:eval_params.nbEval
    out_table(i,:) = out(i).val;
    out_height_vec = [out_height_vec; out(i).data.spike_heights];
    out_burst_freq_vec = [out_burst_freq_vec; 1./out(i).data.bursts_cycle];
    out_spikeburst_vec = [out_spikeburst_vec; out(i).data.bursts_nb_spikes];
end

save(save_path, "loss_prob","out_table", "out_height_vec", "out_burst_freq_vec", "out_spikeburst_vec", "params", "eval_params", "model");


%% Event Loss Plots + Save to CSV

save_path = "stg_loss_neuron1_out.mat";

load(save_path)

window_size = [500 1000];
f = figure;
f.Position = [50 50 window_size];

axs(1) = nexttile;
histogram(out_spikeburst_vec)
ylabel("nb spikes", 'Interpreter', 'tex')

axs(2) = nexttile;
histogram(out_burst_freq_vec)
ylabel("burst freq", 'Interpreter', 'tex')

axs(3) = nexttile;
stem(out_table(:,7), ones(eval_params.nbEval,1))
ylabel("burst length", 'Interpreter', 'tex')

axs(4) = nexttile;
histogram(out_height_vec)
ylabel("spike height", 'Interpreter', 'tex')


T = array2table([out_table(:,3), out_table(:,6), out_table(:,7), out_table(:,9), ones(eval_params.nbEval,1)], 'VariableNames', ["nbspikes", "freq", "length", "height", "val"]);
writetable(T, "stg_loss_avg.csv");

T = array2table([out_spikeburst_vec, ones(length(out_spikeburst_vec),1)], 'VariableNames', ["nbspikes", "val"]);
writetable(T, "stg_loss_nbspike.csv");
T = array2table([out_burst_freq_vec, ones(length(out_burst_freq_vec),1)], 'VariableNames', ["freq", "val"]);
writetable(T, "stg_loss_burstfreq.csv");
T = array2table([out_height_vec, ones(length(out_height_vec),1)], 'VariableNames', ["height", "val"]);
writetable(T, "stg_loss_spikeheight.csv");

%% Event Loss Single Simulation to CSV

model = "simplified_STG";

load_system(model)

t_span = [10, 15];

loss_prob = 0.2;

simin = Simulink.SimulationInput(model);
simin = simin.setVariable('loss_prob', loss_prob, "Workspace", model);
simin = simin.setVariable('a',0.01, "Workspace", model);
out = sim(simin);

values = squeeze(out.events.signals.values);
[red_one, red_time_one] = reduce_points(values(:,1), out.events.time, 0.001);
[red_two, red_time_two] = reduce_points(values(:,2), out.events.time, 0.001);
[red_three, red_time_three] = reduce_points(values(:,3), out.events.time, 0.001);

i_start = find(red_time_one> t_span(1),1)-1;
i_end = find(red_time_one> t_span(2),1);
red_time_one = red_time_one(i_start:i_end);
red_one = red_one(i_start:i_end);
red_spike_time_one = spike_times_finder(red_one, red_time_one);
red_spike_one = ones(size(red_spike_time_one));

i_start = find(red_time_two> t_span(1),1)-1;
i_end = find(red_time_two> t_span(2),1);
red_time_two = red_time_two(i_start:i_end);
red_two = red_two(i_start:i_end);
red_spike_time_two = spike_times_finder(red_two, red_time_two);
red_spike_two = ones(size(red_spike_time_two));

i_start = find(red_time_three> t_span(1),1)-1;
i_end = find(red_time_three> t_span(2),1);
red_time_three = red_time_three(i_start:i_end);
red_three = red_three(i_start:i_end);
red_spike_time_three = spike_times_finder(red_three, red_time_three);
red_spike_three = ones(size(red_spike_time_three));

figure;
hold on
stem(red_spike_time_one, red_spike_one, "Color",ulgColors.BlueGreen.hex);
stem(red_spike_time_two, red_spike_two, "Color",ulgColors.YellowOrange.hex);
stem(red_spike_time_three, red_spike_three, "Color",ulgColors.Red.hex);
xlim(t_span)

red_spike_time_one = red_spike_time_one - t_span(1);
red_spike_time_two = red_spike_time_two - t_span(1);
red_spike_time_three = red_spike_time_three - t_span(1);

T = array2table([red_spike_time_one, red_spike_one], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron1_event_loss.csv");
T = array2table([red_spike_time_two, red_spike_two], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron2_event_loss.csv");
T = array2table([red_spike_time_three, red_spike_three], 'VariableNames', ["time", "value"]);
writetable(T, "STG_neuron3_event_loss.csv");

%% Script Local Functions

function spike_times = spike_times_finder(spike_trace, spike_trace_t, tresh)
    if nargin < 3
        tresh = 0.5;
    end
    spike_trace_length = length(spike_trace);
    spike_times = zeros(spike_trace_length/2 + 1,1);
    spike_ind = 0;
    for i = 2:spike_trace_length
        if spike_trace(i-1) <= tresh && spike_trace(i) > tresh
            spike_ind = spike_ind + 1;
            spike_times(spike_ind) = (spike_trace_t(i-1) + spike_trace_t(i))/2;
        end
    end
    spike_times = spike_times(1:spike_ind);
end

function [X, T] = clip_to_span(X, T, span)
    N = length(T); 
    if span(2) < span(1)
        error("Span must be increasing")
    end
    if T(1) > span(1)
        error("Start of data must be before span")
    end
    if T(N) < span(2)
        error("End of data must be after span")
    end
 
    X(1) = X(1) * (T(2) - span(1))/(T(2) - T(1)) + X(2) * (span(1) - T(1))/(T(2) - T(1));
    T(1) = span(1);

    X(N) = X(N-1) * (T(N) - span(2))/(T(N) - T(N-1)) + X(N) * (span(2) - T(N-1))/(T(N) - T(N-1));
    T(N) = span(2);
end