model = "synapse_example_circuit";

load_system(model)

t_span = [11.3,12.9];

out = sim(model);

values = squeeze(out.synapses.signals.values);
[red_one, red_time_one] = reduce_points(values(:,1), out.synapses.time, 0.001);
[red_two, red_time_two] = reduce_points(values(:,2), out.synapses.time, 0.001);
[red_three, red_time_three] = reduce_points(values(:,3), out.synapses.time, 0.001);
[red_four, red_time_four] = reduce_points(values(:,4), out.synapses.time, 0.001);

i_start = find(red_time_one> t_span(1),1)-1;
i_end = find(red_time_one> t_span(2),1);
red_time_one = red_time_one(i_start:i_end);
red_one = red_one(i_start:i_end);

i_start = find(red_time_two> t_span(1),1)-1;
i_end = find(red_time_two> t_span(2),1);
red_time_two = red_time_two(i_start:i_end);
red_two = red_two(i_start:i_end);

i_start = find(red_time_three> t_span(1),1)-1;
i_end = find(red_time_three> t_span(2),1);
red_time_three = red_time_three(i_start:i_end);
red_three = red_three(i_start:i_end);

i_start = find(red_time_four> t_span(1),1)-1;
i_end = find(red_time_four> t_span(2),1);
red_time_four = red_time_four(i_start:i_end);
red_four = red_four(i_start:i_end);

[red_one, red_time_one] = clip_to_span(red_one, red_time_one, t_span);
[red_two, red_time_two] = clip_to_span(red_two, red_time_two, t_span);
[red_three, red_time_three] = clip_to_span(red_three, red_time_three, t_span);
[red_four, red_time_four] = clip_to_span(red_four, red_time_four, t_span);

red_four_spike_time = spike_times_finder(red_four, red_time_four);
red_four_spike_value = ones(size(red_four_spike_time));


figure;
hold on
plot(red_time_one, red_one, "Color",ulgColors.BlueGreen.hex);
plot(red_time_two, red_two, "Color",ulgColors.YellowOrange.hex);
plot(red_time_three, red_three, "Color",ulgColors.Red.hex);
stem(red_four_spike_time, red_four_spike_value, "Color",ulgColors.Violet.hex);
xlim(t_span)

red_time_one = red_time_one - t_span(1);
red_time_two = red_time_two - t_span(1);
red_time_three = red_time_three - t_span(1);
red_four_spike_time = red_four_spike_time - t_span(1);

T = array2table([red_time_one, red_one], 'VariableNames', ["time", "value"]);
writetable(T, "example_synapse1.csv");
T = array2table([red_time_two, red_two], 'VariableNames', ["time", "value"]);
writetable(T, "example_synapse2.csv");
T = array2table([red_time_three, red_three], 'VariableNames', ["time", "value"]);
writetable(T, "example_synapse3.csv");
T = array2table([red_four_spike_time, red_four_spike_value], 'VariableNames', ["time", "value"]);
writetable(T, "example_synapse4_spike.csv");


%% 

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