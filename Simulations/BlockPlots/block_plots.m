%% Bistable Relay + Threshold Detector Plots

model = "bistable_and_rebound_plots";

load_system(model)

out = sim(model);

t_span = [1, 3.75];

values = squeeze(out.bistable.signals.values);
[red_low, red_time_low] = reduce_points(values(:,1), out.bistable.time, 0.001);
[red_tanh, red_time_tanh] = reduce_points(values(:,2), out.bistable.time, 0.001);
[red_in_pos, red_time_in_pos] = reduce_points(values(:,3), out.bistable.time, 0.001);
[red_in_neg, red_time_in_neg] = reduce_points(values(:,4), out.bistable.time, 0.001);

i_start = find(red_time_low> t_span(1),1)-1;
i_end = find(red_time_low> t_span(2),1);
red_time_low = red_time_low(i_start:i_end);
red_low = red_low(i_start:i_end);

i_start = find(red_time_tanh> t_span(1),1)-1;
i_end = find(red_time_tanh> t_span(2),1);
red_time_tanh = red_time_tanh(i_start:i_end);
red_tanh = red_tanh(i_start:i_end);

i_start = find(red_time_in_pos> t_span(1),1)-1;
i_end = find(red_time_in_pos> t_span(2),1);
red_time_in_pos = red_time_in_pos(i_start:i_end);
red_in_pos = red_in_pos(i_start:i_end);

i_start = find(red_time_in_neg> t_span(1),1)-1;
i_end = find(red_time_in_neg> t_span(2),1);
red_time_in_neg = red_time_in_neg(i_start:i_end);
red_in_neg = red_in_neg(i_start:i_end);

red_in_pos_spike_time = spike_times_finder(red_in_pos, red_time_in_pos);
red_in_pos_spike_value = ones(size(red_in_pos_spike_time));
red_in_neg_spike_time = spike_times_finder(red_in_neg, red_time_in_neg);
red_in_neg_spike_value = ones(size(red_in_neg_spike_time));

[red_low, red_time_low] = clip_to_span(red_low, red_time_low, t_span);
[red_tanh, red_time_tanh] = clip_to_span(red_tanh, red_time_tanh, t_span);
[red_in_pos, red_time_in_pos] = clip_to_span(red_in_pos, red_time_in_pos, t_span);
[red_in_neg, red_time_in_neg] = clip_to_span(red_in_neg, red_time_in_neg, t_span);

figure;
hold on
plot(red_time_low, red_low, "Color",ulgColors.BlueGreen.hex);
plot(red_time_tanh, red_tanh, "Color",ulgColors.Red.hex);
plot(red_time_in_pos, red_in_pos, "Color",ulgColors.Violet.hex);
stem(red_in_pos_spike_time, red_in_pos_spike_value, "Color",ulgColors.Yellow.hex);
plot(red_time_in_neg, -red_in_neg, "Color",ulgColors.Lavander.hex);
stem(red_in_neg_spike_time, -red_in_neg_spike_value, "Color",ulgColors.YellowOrange.hex);
xlim(t_span)

red_time_low = red_time_low - t_span(1);
red_time_tanh = red_time_tanh - t_span(1);
red_time_in_pos = red_time_in_pos - t_span(1);
red_in_pos_spike_time = red_in_pos_spike_time - t_span(1);
red_time_in_neg = red_time_in_neg - t_span(1);
red_in_neg_spike_time = red_in_neg_spike_time - t_span(1);

T = array2table([red_time_low, red_low], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_lowpass.csv");
T = array2table([red_time_tanh, red_tanh], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_tanh.csv");
T = array2table([red_time_in_pos, red_in_pos], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_input_pos.csv");
T = array2table([red_in_pos_spike_time, red_in_pos_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_input_pos_spike.csv");
T = array2table([red_time_in_neg, red_in_neg], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_input_neg.csv");
T = array2table([red_in_neg_spike_time, red_in_neg_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "bistable_input_neg_spike.csv");


t_span = [0.75, 1.75];

values = squeeze(out.rebound.signals.values);
[red_out, red_time_out] = reduce_points(values(:,1), out.rebound.time, 0.0003);
[red_inhib, red_time_inhib] = reduce_points(values(:,2), out.rebound.time, 0.0003);
[red_in, red_time_in] = reduce_points(values(:,3), out.rebound.time, 0.0003);

i_start = find(red_time_out> t_span(1),1)-1;
i_end = find(red_time_out> t_span(2),1);
red_time_out = red_time_out(i_start:i_end);
red_out = red_out(i_start:i_end);

i_start = find(red_time_inhib> t_span(1),1)-1;
i_end = find(red_time_inhib> t_span(2),1);
red_time_inhib = red_time_inhib(i_start:i_end);
red_inhib = red_inhib(i_start:i_end);

i_start = find(red_time_in> t_span(1),1)-1;
i_end = find(red_time_in> t_span(2),1);
red_time_in = red_time_in(i_start:i_end);
red_in = red_in(i_start:i_end);

red_out_spike_time = spike_times_finder(red_out, red_time_out);
red_out_spike_value = ones(size(red_out_spike_time));

[red_out, red_time_out] = clip_to_span(red_out, red_time_out, t_span);
[red_inhib, red_time_inhib] = clip_to_span(red_inhib, red_time_inhib, t_span);
[red_in, red_time_in] = clip_to_span(red_in, red_time_in, t_span);

figure;
hold on
plot(red_time_out, red_out, "Color", ulgColors.BlueGreen.hex);
stem(red_out_spike_time, red_out_spike_value, "Color",ulgColors.Yellow.hex);
plot(red_time_inhib, red_inhib, "Color",ulgColors.Red.hex);
plot(red_time_in, red_in, "Color",ulgColors.Violet.hex);
xlim(t_span)

red_time_out = red_time_out - t_span(1);
red_time_inhib = red_time_inhib - t_span(1);
red_time_in = red_time_in - t_span(1);
red_out_spike_time = red_out_spike_time - t_span(1);

T = array2table([red_time_out, red_out], 'VariableNames', ["time", "value"]);
writetable(T, "rebound_output.csv");
T = array2table([red_out_spike_time, red_out_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "rebound_output_spike.csv");
T = array2table([red_time_inhib, red_inhib], 'VariableNames', ["time", "value"]);
writetable(T, "rebound_inhib.csv");
T = array2table([red_time_in, red_in], 'VariableNames', ["time", "value"]);
writetable(T, "rebound_input.csv");

%% Neuron Plot
 
model = "neuron_plots";

load_system(model)

out = sim(model);

t_span = [4.9, 5.35];

values = squeeze(out.neuron_behavior.signals.values);
[red_input, red_time_input] = reduce_points(values(:,3), out.neuron_behavior.time, 0.001);
[red_event, red_time_event] = reduce_points(values(:,1), out.neuron_behavior.time, 0.001);
[red_output, red_time_output] = reduce_points(values(:,2), out.neuron_behavior.time, 0.001);

i_start = find(red_time_input> t_span(1),1)-1;
i_end = find(red_time_input> t_span(2),1);
red_time_input = red_time_input(i_start:i_end);
red_input = red_input(i_start:i_end);

i_start = find(red_time_event> t_span(1),1)-1;
i_end = find(red_time_event> t_span(2),1);
red_time_event = red_time_event(i_start:i_end);
red_event = red_event(i_start:i_end);

i_start = find(red_time_output> t_span(1),1)-1;
i_end = find(red_time_output> t_span(2),1);
red_time_output = red_time_output(i_start:i_end);
red_output = red_output(i_start:i_end);

red_event_spike_time = spike_times_finder(red_event, red_time_event);
red_event_spike_value = ones(size(red_event_spike_time));

[red_input, red_time_input] = clip_to_span(red_input, red_time_input, t_span);
[red_event, red_time_event] = clip_to_span(red_output, red_time_output, t_span);
[red_output, red_time_output] = clip_to_span(red_output, red_time_output, t_span);

figure;
hold on
plot(red_time_input, red_input, "Color",ulgColors.BlueGreen.hex);
plot(red_time_event, red_event, "Color",ulgColors.YellowOrange.hex);
plot(red_time_output, red_output, "Color",ulgColors.Red.hex);
stem(red_event_spike_time, red_event_spike_value, "Color",ulgColors.Violet.hex);
xlim(t_span)

red_time_input = red_time_input - t_span(1);
red_time_output = red_time_output - t_span(1);
red_time_event = red_time_event - t_span(1);
red_event_spike_time = red_event_spike_time - t_span(1);

T = array2table([red_time_input, red_input], 'VariableNames', ["time", "value"]);
writetable(T, "neuron_input.csv");
T = array2table([red_time_event, red_event], 'VariableNames', ["time", "value"]);
writetable(T, "neuron_output_event.csv");
T = array2table([red_event_spike_time, red_event_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "neuron_output_spike.csv");
T = array2table([red_time_output, red_output], 'VariableNames', ["time", "value"]);
writetable(T, "neuron_output.csv");

%% Additive Synapses Plots

model = "synapse_plots";

load_system(model)

out = sim(model);

t_span = [3.35, 3.95];

values = squeeze(out.depressing.signals.values);
[red_input, red_time_input] = reduce_points(values(:,1), out.depressing.time, 0.001);
[red_output, red_time_output] = reduce_points(values(:,2), out.depressing.time, 0.001);

i_start = find(red_time_input> t_span(1),1)-1;
i_end = find(red_time_input> t_span(2),1);
red_time_input = red_time_input(i_start:i_end);
red_input = red_input(i_start:i_end);

i_start = find(red_time_output> t_span(1),1)-1;
i_end = find(red_time_output> t_span(2),1);
red_time_output = red_time_output(i_start:i_end);
red_output = red_output(i_start:i_end);

red_input_spike_time = spike_times_finder(red_input, red_time_input);
red_input_spike_value = ones(size(red_input_spike_time));

[red_input, red_time_input] = clip_to_span(red_input, red_time_input, t_span);
[red_output, red_time_output] = clip_to_span(red_output, red_time_output, t_span);

figure;
hold on
plot(red_time_input, red_input, "Color",ulgColors.BlueGreen.hex);
plot(red_time_output, red_output, "Color",ulgColors.Red.hex);
stem(red_input_spike_time, red_input_spike_value, "Color",ulgColors.Violet.hex);
xlim(t_span)

red_time_input = red_time_input - t_span(1);
red_time_output = red_time_output - t_span(1);
red_input_spike_time = red_input_spike_time - t_span(1);

T = array2table([red_time_input, red_input], 'VariableNames', ["time", "value"]);
writetable(T, "depressing_input.csv");
T = array2table([red_input_spike_time, red_input_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "depressing_input_spike.csv");
T = array2table([red_time_output, red_output], 'VariableNames', ["time", "value"]);
writetable(T, "depressing_output.csv");


values = squeeze(out.facilitating.signals.values);
[red_input, red_time_input] = reduce_points(values(:,1), out.facilitating.time, 0.001);
[red_output, red_time_output] = reduce_points(values(:,2), out.facilitating.time, 0.001);

i_start = find(red_time_input> t_span(1),1)-1;
i_end = find(red_time_input> t_span(2),1);
red_time_input = red_time_input(i_start:i_end);
red_input = red_input(i_start:i_end);

i_start = find(red_time_output> t_span(1),1)-1;
i_end = find(red_time_output> t_span(2),1);
red_time_output = red_time_output(i_start:i_end);
red_output = red_output(i_start:i_end);

red_input_spike_time = spike_times_finder(red_input, red_time_input);
red_input_spike_value = ones(size(red_input_spike_time));

[red_input, red_time_input] = clip_to_span(red_input, red_time_input, t_span);
[red_output, red_time_output] = clip_to_span(red_output, red_time_output, t_span);

figure;
hold on
plot(red_time_input, red_input, "Color",ulgColors.BlueGreen.hex);
plot(red_time_output, red_output, "Color",ulgColors.Red.hex);
stem(red_input_spike_time, red_input_spike_value, "Color",ulgColors.Violet.hex);
xlim(t_span)

red_time_input = red_time_input - t_span(1);
red_time_output = red_time_output - t_span(1);
red_input_spike_time = red_input_spike_time - t_span(1);

T = array2table([red_time_input, red_input], 'VariableNames', ["time", "value"]);
writetable(T, "facilitating_input.csv");
T = array2table([red_input_spike_time, red_input_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "facilitating_input_spike.csv");
T = array2table([red_time_output, red_output], 'VariableNames', ["time", "value"]);
writetable(T, "facilitating_output.csv");

%% Modulation Plots

model = "modulation_plots";

load_system(model)

out = sim(model);

t_span = [0, 10];

values = squeeze(out.modulation.signals.values);
[red_output, red_time_output] = reduce_points(values(:,1), out.modulation.time, 0.001);
[red_input_pos, red_time_input_pos] = reduce_points(values(:,2), out.modulation.time, 0.001);
[red_input_neg, red_time_input_neg] = reduce_points(values(:,3), out.modulation.time, 0.001);


i_start = find(red_time_output> t_span(1),1)-1;
i_end = find(red_time_output> t_span(2),1);
red_time_output = red_time_output(i_start:i_end);
red_output = red_output(i_start:i_end);

i_start = find(red_time_input_pos> t_span(1),1)-1;
i_end = find(red_time_input_pos> t_span(2),1);
red_time_input_pos = red_time_input_pos(i_start:i_end);
red_input_pos = red_input_pos(i_start:i_end);

i_start = find(red_time_input_neg> t_span(1),1)-1;
i_end = find(red_time_input_neg> t_span(2),1);
red_time_input_neg = red_time_input_neg(i_start:i_end);
red_input_neg = red_input_neg(i_start:i_end);



red_input_pos_spike_time = spike_times_finder(red_input_pos, red_time_input_pos);
red_input_pos_spike_value = ones(size(red_input_pos_spike_time));

red_input_neg_spike_time = spike_times_finder(red_input_neg, red_time_input_neg);
red_input_neg_spike_value = ones(size(red_input_neg_spike_time));

[red_output, red_time_output] = clip_to_span(red_output, red_time_output, t_span);
[red_input_pos, red_time_input_pos] = clip_to_span(red_input_pos, red_time_input_pos, t_span);
[red_input_neg, red_time_input_neg] = clip_to_span(red_input_neg, red_time_input_neg, t_span);

figure;
hold on
plot(red_time_input_pos, red_input_pos, "Color",ulgColors.BlueGreen.hex);
stem(red_input_pos_spike_time, red_input_pos_spike_value, "Color",ulgColors.Yellow.hex);
plot(red_time_input_neg, red_input_neg, "Color",ulgColors.Violet.hex);
stem(red_input_neg_spike_time, -red_input_neg_spike_value, "Color",ulgColors.Yellow.hex);
plot(red_time_output, red_output, "Color",ulgColors.Red.hex);
xlim(t_span)

red_time_input_pos = red_time_input_pos - t_span(1);
red_time_input_neg = red_time_input_neg - t_span(1);
red_time_output = red_time_output - t_span(1);

T = array2table([red_time_input_pos, red_input_pos], 'VariableNames', ["time", "value"]);
writetable(T, "modulation_input_pos.csv");
T = array2table([red_input_pos_spike_time, red_input_pos_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "modulation_input_pos_spike.csv");
T = array2table([red_time_input_neg, red_input_neg], 'VariableNames', ["time", "value"]);
writetable(T, "modulation_input_neg.csv");
T = array2table([red_input_neg_spike_time, red_input_neg_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "modulation_input_neg_spike.csv");
T = array2table([red_time_output, red_output], 'VariableNames', ["time", "value"]);
writetable(T, "modulation_output.csv");

%% Event-Analog Signal Multiplier Plots

model = "eventdemux_plots";

load_system(model)

out = sim(model);


values = squeeze(out.eventdemux.signals.values);
[red_output_pos, red_time_output_pos] = reduce_points(values(:,1), out.eventdemux.time, 0.001);
[red_output_neg, red_time_output_neg] = reduce_points(values(:,2), out.eventdemux.time, 0.001);
[red_input, red_time_input] = reduce_points(values(:,3), out.eventdemux.time, 0.001);
[red_signal, red_time_signal] = reduce_points(values(:,4), out.eventdemux.time, 0.001);


red_output_pos_spike_time = spike_times_finder(red_output_pos, red_time_output_pos);
red_output_pos_spike_value = ones(size(red_output_pos_spike_time));
red_output_neg_spike_time = spike_times_finder(red_output_neg, red_time_output_neg);
red_output_neg_spike_value = ones(size(red_output_neg_spike_time));
red_input_spike_time = spike_times_finder(red_input, red_time_input);
red_input_spike_value = ones(size(red_input_spike_time));

figure;
hold on
plot(red_time_input, red_input, "Color",ulgColors.BlueGreen.hex);
plot(red_time_signal, red_signal, "Color",ulgColors.Violet.hex);
stem(red_input_spike_time, 0.2*red_input_spike_value, "Color",ulgColors.Lavander.hex);
plot(red_time_output_pos, red_output_pos, "Color",ulgColors.Red.hex);
stem(red_output_pos_spike_time, red_output_pos_spike_value, "Color",ulgColors.Yellow.hex);
plot(red_time_output_neg, -red_output_neg, "Color",ulgColors.Orange.hex);
stem(red_output_neg_spike_time, -red_output_neg_spike_value, "Color",ulgColors.YellowOrange.hex);
%xlim(t_span)


T = array2table([red_time_output_pos, red_output_pos], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_output_pos.csv");
T = array2table([red_output_pos_spike_time, red_output_pos_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_output_pos_spike.csv");
T = array2table([red_time_output_neg, red_output_neg], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_output_neg.csv");
T = array2table([red_output_neg_spike_time, red_output_neg_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_output_neg_spike.csv");
T = array2table([red_time_input, red_input], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_input.csv");
T = array2table([red_input_spike_time, red_input_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_input_spike.csv");
T = array2table([red_time_signal, red_signal], 'VariableNames', ["time", "value"]);
writetable(T, "eventdemux_signal.csv");


%% 

function spike_times = spike_times_finder(spike_trace, spike_trace_t, tresh)
    if nargin < 3
        tresh = 0.5;
    end
    spike_trace_length = length(spike_trace);
    spike_times = zeros(ceil(spike_trace_length/2) + 1,1);
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
