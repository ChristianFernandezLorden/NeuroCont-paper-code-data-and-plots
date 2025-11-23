function newout = postsim_neuronal_signal_analysis(o, eval_params, params)
    newout = struct();
    
    newout.data.spike_heights = [];
    newout.data.bursts_cycle = [];
    newout.data.bursts_size = [];
    newout.data.bursts_nb_spikes = [];
    newout.data.duty = [];
    newout.data.spike_cycle = [];

    try 
        startAnalyseTime = eval_params.StartAnalyseTime;
        
        useful = find(o.voltage.time >= startAnalyseTime);
        T = o.voltage.time(useful);
        V = o.voltage.signals.values(useful, 1);
    
        nb_spikes = NaN;
        b_length = NaN;
        duty = NaN;
        p_length = NaN;
        intra_freq = NaN;
        inter_freq = NaN;
        spike_height = NaN;
    
        prev = V(1);
        up_crossings = zeros(length(T), 1);
        down_crossings = zeros(length(T), 1);
        ind_upcross = 0;
        ind_downcross = 0;
    
        for i = 2:length(T)
            curr = V(i);
            if prev < 0 && curr > 0
                ind_upcross = ind_upcross + 1;
                % Linear interpolation
                t_prev = T(i-1);
                t_curr = T(i);
                t_inter = (t_prev*curr - t_curr*prev) / (curr - prev);
                up_crossings(ind_upcross) = t_inter;
            elseif prev > 0 && curr < 0
                ind_downcross = ind_downcross + 1;
                % Linear interpolation
                t_prev = T(i-1);
                t_curr = T(i);
                t_inter = (t_prev*curr - t_curr*prev) / (curr - prev);
                down_crossings(ind_downcross) = t_inter;
            end
            prev = curr;
        end
        up_crossings = up_crossings(1:ind_upcross);
        down_crossings = down_crossings(1:ind_downcross);
        
        if length(up_crossings) > 3 && length(down_crossings) > 3
            % down separation -> between a down and an up crossing
            % up separation -> between an up and a down crossing

            % Alternance up -> down -> up -> down ... (up -> down) or (up-> down -> up)
            if up_crossings(1) < down_crossings(1)
                down_sep = up_crossings(2:end) - down_crossings(1:length(up_crossings)-1);
                up_sep = down_crossings(1:end) - up_crossings(1:length(down_crossings));
                
                spike_heights = zeros(length(down_crossings), 1);
                for i = 1:length(down_crossings)
                    ind_start =  find(T>up_crossings(i),1);
                    if isempty(ind_start)
                        ind_start = 1;
                    end
                    ind_end = find(T<down_crossings(i),1, 'last');
                    if isempty(ind_end)
                        ind_end = length(V);
                    end
                    spike_heights(i) = max(V(ind_start:ind_end)) - min(V(ind_start:ind_end));
                end
            % Alternance down -> up -> down ... (down-> up -> down) or (down -> up)
            else
                down_sep = up_crossings(1:end) - down_crossings(1:length(up_crossings));
                up_sep = down_crossings(2:end) - up_crossings(1:length(down_crossings)-1);

                spike_heights = zeros(length(up_crossings), 1);
                for i = 1:length(down_crossings)-1
                    ind_start =  find(T>up_crossings(i),1);
                    if isempty(ind_start)
                        ind_start = 1;
                    end
                    ind_end = find(T<down_crossings(i+1),1, 'last');
                    if isempty(ind_end)
                        ind_end = length(V);
                    end
                    spike_heights(i) = max(V(ind_start:ind_end)) - min(V(ind_start:ind_end));
                end
            end
            newout.data.spike_heights = spike_heights;

            spike_height = mean(spike_heights);
            %cross_sep = crossings(2:end) - crossings(1:end-1);
            
            % Test bursting on the down separation (tile between burst is
            % negative)
            [down_sep_sort, idx_down] = sort(down_sep);
            ind_min_down = 1;
            num_min_down = 1;
            sum_min_down = down_sep_sort(1);
            ind_max_down = length(down_sep_sort);
            sum_max_down = down_sep_sort(end);
            num_max_down = 1;
            while ind_min_down < ind_max_down-1
                diff_min = down_sep_sort(ind_min_down+1) - sum_min_down/num_min_down;
                diff_max = sum_max_down/num_max_down - down_sep_sort(ind_max_down-1);
                if diff_min < diff_max
                    ind_min_down = ind_min_down + 1;
                    sum_min_down = sum_min_down + down_sep_sort(ind_min_down);
                    num_min_down = num_min_down + 1;
                else
                    ind_max_down = ind_max_down - 1;
                    sum_max_down = sum_max_down + down_sep_sort(ind_max_down);
                    num_max_down = num_max_down + 1;
                end
            end
    
            [up_sep_sort, idx_up] = sort(up_sep);
            ind_min_up = 1;
            num_min_up = 1;
            sum_min_up = up_sep_sort(1);
            ind_max_up = length(up_sep_sort);
            sum_max_up = up_sep_sort(end);
            num_max_up = 1;
            while ind_min_up < ind_max_up-1
                diff_min = up_sep_sort(ind_min_up+1) - sum_min_up/num_min_up;
                diff_max = sum_max_up/num_max_up - up_sep_sort(ind_max_up-1);
                if diff_min < diff_max
                    ind_min_up = ind_min_up + 1;
                    sum_min_up = sum_min_up + up_sep_sort(ind_min_up);
                    num_min_up = num_min_up + 1;
                else
                    ind_max_up = ind_max_up - 1;
                    sum_max_up = sum_max_up + up_sep_sort(ind_max_up);
                    num_max_up = num_max_up + 1;
                end
            end
    
            % Criterion of bursting
            if sum_max_down/num_max_down > 4*sum_min_down/num_min_down
                burst_gaps = zeros(size(down_sep));
                burst_gaps(idx_down(ind_max_down:end)) = 1;
                burst_gaps = logical(burst_gaps);
                
                % Change scope of V for power computation
                first_gap = find(burst_gaps, 1, 'first');
                % up then down
                if up_crossings(1) < down_crossings(1)
                    T_start = down_crossings(1) + sum(up_sep(2:first_gap)) + sum(down_sep(1:first_gap-1));
                else
                    T_start = down_crossings(1) + sum(up_sep(1:first_gap-1)) + sum(down_sep(1:first_gap-1));
                end
                last_gap = find(burst_gaps, 1, 'last');
                if up_crossings(1) < down_crossings(1)
                    T_end = down_crossings(1) + sum(up_sep(2:last_gap)) + sum(down_sep(1:last_gap-1));
                else
                    T_end = down_crossings(1) + sum(up_sep(1:last_gap-1)) + sum(down_sep(1:last_gap-1));
                end
                id_start = find(T <= T_start, 1, 'last');
                id_end = find(T <= T_end, 1, 'last') + 1;
                V = V(id_start:id_end);
                
                % Same condition as bursting for plateau
                if sum_max_up/num_max_up > 4*sum_min_up/num_min_up
                    type = "plateau_bursting";
                    plateau_gaps = zeros(size(up_sep));
                    plateau_gaps(idx_up(ind_max_up:end)) = 1;
                    plateau_gaps = logical(plateau_gaps);
                    p_length = mean(up_sep(plateau_gaps));
                else
                    type = "bursting";
                end
    
                bursts_cycle = zeros(sum(burst_gaps, 'all'), 1);
                bursts_size = zeros(sum(burst_gaps, 'all'), 1);
                bursts_nb_spikes = zeros(sum(burst_gaps, 'all'), 1);
                bursts_duty = zeros(sum(burst_gaps, 'all'), 1);
                spike_cycle = zeros(sum(burst_gaps, 'all'), 1);
                
                if up_crossings(1) < down_crossings(1)
                    up_sep = up_sep(2:end);
                end
    
                min_length = min(length(up_sep), length(down_sep));
                down_sep = down_sep(1:min_length);
                up_sep = up_sep(1:min_length);
    
                ind_start = 0;
                for i = 1:min_length
                    if burst_gaps(i)
                        ind_start = i+1;
                        break;
                    end
                end
                ind_burst = 0;
                cumulative_time = up_sep(ind_start-1);
                cumulative_spike_cycle = 0;
                nb_spike = 0;
                for i = ind_start:min_length
                    nb_spike = nb_spike+1;
                    % When a new bursting gap is found the cycle is done
                    if burst_gaps(i)
                        ind_burst = ind_burst+1;
    
                        bursts_nb_spikes(ind_burst) = nb_spike;
                        bursts_cycle(ind_burst) = cumulative_time + down_sep(i);
                        bursts_size(ind_burst) = cumulative_time;
                        bursts_duty(ind_burst) = bursts_size(ind_burst)/bursts_cycle(ind_burst);
                        spike_cycle(ind_burst) = cumulative_spike_cycle/(nb_spike-1);
    
                        cumulative_time = up_sep(i);
                        cumulative_spike_cycle = 0;
                        nb_spike = 0;
                    else
                        cumulative_time = cumulative_time + down_sep(i) + up_sep(i);
                        cumulative_spike_cycle = cumulative_spike_cycle + down_sep(i) + up_sep(i);
                    end
                end
    
                bursts_cycle = bursts_cycle(1:ind_burst);
                bursts_size = bursts_size(1:ind_burst);
                bursts_nb_spikes = bursts_nb_spikes(1:ind_burst);
                bursts_duty = bursts_duty(1:ind_burst);
                spike_cycle = spike_cycle(1:ind_burst);

                newout.data.bursts_cycle = bursts_cycle;
                newout.data.bursts_size = bursts_size;
                newout.data.bursts_nb_spikes = bursts_nb_spikes;
                newout.data.duty = bursts_duty;
                newout.data.spike_cycle = spike_cycle;
    
                nb_spikes = mean(bursts_nb_spikes);
                b_length = mean(bursts_size);
                intra_freq = 1/mean(spike_cycle);
                inter_freq = 1/mean(bursts_cycle);
                duty = mean(bursts_duty);
                
                
            % If not bursting then spiking
            else
                type = "spiking";
                nb_spikes = 1;
                
                % Up then down (more logical)
                if up_crossings(1) > down_crossings(1)
                    down_sep = down_sep(2:end);
                end
                
                % Change scope of V for power computation
                T_start = up_crossings(1);
                T_end = up_crossings(end);
                id_start = find(T <= T_start, 1, 'last');
                id_end = find(T <= T_end, 1, 'last') + 1;
                V = V(id_start:id_end);
                
                min_length = min(length(up_sep), length(down_sep));
                down_sep = down_sep(1:min_length);
                up_sep = up_sep(1:min_length);
                total_sep = down_sep + up_sep;

                newout.data.spike_cycle = total_sep;
                newout.data.duty = up_sep./total_sep;
                
                intra_freq = 1/mean(total_sep);
                inter_freq = NaN;
                duty = mean(up_sep./total_sep);
            end
    
        else
    
            if mean(V) < 0
                type = "silent_hyperpolarized";
            else
                type = "silent_depolarized";
            end
        end
    
        power_signal = V;
        power_signal(power_signal < 0.0) = 0.0;
        power = mean(abs(power_signal));
        
        type_dict = containers.Map( ... % Order of power
            ["silent_hyperpolarized", "spiking", "bursting", "plateau_bursting", "silent_depolarized"], ...
            [0, 1, 2, 3, 4]...
        );
        %newout = struct("val", [type_dict(type), power, nb_spikes, duty, intra_freq, inter_freq, b_length, p_length, spike_height]);
        newout.val = [type_dict(type), power, nb_spikes, duty, intra_freq, inter_freq, b_length, p_length, spike_height];
        return
    catch E
        %newout = struct("val", [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]);
        newout.val = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        return
    end
end