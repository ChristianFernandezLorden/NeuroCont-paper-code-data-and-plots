model = 'neuron_exc_type';

load_system(model)

signals = {'typeI', 'typeII', 'typeIII', 'typeB'};
tspans = {[1.75, 4.75], [1.75, 4.75], [1.75, 3.25], [1.75, 4.5]};

out = sim(model);

for i = 1:length(signals)
    signal_name = signals{i};
    tspan = tspans{i};
    
    [output, output_T] = reduce_points(out.(signal_name).signals.values(:,1), out.(signal_name).time, 0.001);
    [output, output_T] = clip_to_span(output, output_T, tspan);
    [input, input_T] = reduce_points(out.(signal_name).signals.values(:,2), out.(signal_name).time, 0.001);
    [input, input_T] = clip_to_span(input, input_T, tspan);
    
    f = figure;
    t = tiledlayout(3,1);
    nexttile([2,1])
    plot(output_T, output);
    nexttile
    plot(input_T, input);
    
    
    T = array2table([output_T, output], 'VariableNames', ["time", "value"]);
    writetable(T, ['exc_types_' signal_name '_output.csv']);
    T = array2table([input_T, input], 'VariableNames', ["time", "value"]);
    writetable(T, ['exc_types_' signal_name '_input.csv']);
end


%%

function [X, T] = clip_to_span(X, T, span)
    N = length(T); 
    if span(2) < span(1)
        error("Span must be increasing")
    end
    if T(1) >= span(1)
        error("Start of data must be before span")
    end
    if T(N) <= span(2)
        error("End of data must be after span")
    end


    i_before = find(T > span(1),1)-1;
    i_after = find(T > span(2),1);
    T = T(i_before:i_after);
    X = X(i_before:i_after);
    N = length(T); 
 
    X(1) = X(1) * (T(2) - span(1))/(T(2) - T(1)) + X(2) * (span(1) - T(1))/(T(2) - T(1));
    T(1) = span(1);

    X(N) = X(N-1) * (T(N) - span(2))/(T(N) - T(N-1)) + X(N) * (span(2) - T(N-1))/(T(N) - T(N-1));
    T(N) = span(2);

    T = T - span(1);
end