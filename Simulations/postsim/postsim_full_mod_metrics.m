function newout = postsim_full_mod_metrics(o, eval_params, params)
 % Speed metrcis ?
    newout = struct();
    try  % Try-Catch of reading end params
        startAnalyseTime = eval_params.StartAnalyseTime;

        
        filter_index = find(o.gsyn.time >= startAnalyseTime);
        %T_gsyn = o.gsyn.time(filter_index);
        gsyn = squeeze(o.gsyn.signals.values);
        gsyn_pos = gsyn(filter_index, 1);
        gsyn_neg = gsyn(filter_index, 2);

        [std_gsyn_pos, mean_gsyn_pos] = std(gsyn_pos);
        [std_gsyn_neg, mean_gsyn_neg] = std(gsyn_neg);

        filter_index = find(o.gsm_mot.time >= startAnalyseTime);
        %T_gsm_mot = o.gsm_mot.time(filter_index);
        gsm_mot = squeeze(o.gsm_mot.signals.values);
        gsm_mot_pos = gsm_mot(filter_index, 1);
        gsm_mot_neg = gsm_mot(filter_index, 2);
        
        [std_gsm_mot_pos, mean_gsm_mot_pos] = std(gsm_mot_pos);
        [std_gsm_mot_neg, mean_gsm_mot_neg] = std(gsm_mot_neg);
    catch E
        % Return whats possible (avoid crashes)
        newout.val =   [NaN, NaN, NaN, NaN,...
                        NaN, NaN, NaN, NaN,...
                        NaN, NaN, NaN, NaN,...
                        NaN, NaN, NaN, NaN];
        return
    end
        
        
    try  % Try-Catch of peak detection
        % === Find Mean and Std of Amplitude ===
        
        % Extract Data  
        filter_index = find(o.measure.time >= startAnalyseTime);
        T_angle_mes = o.measure.time(filter_index);
        angle_mes = squeeze(o.measure.signals.values);
        angle_mes = angle_mes(filter_index, 1);
        
        
        [max_amp, max_amp_T] = findpeaks(angle_mes,T_angle_mes, ...
                            'MinPeakProminence', eval_params.peakProminence);
        [min_amp, min_amp_T] = findpeaks(-angle_mes,T_angle_mes, ...
                            'MinPeakProminence', eval_params.peakProminence);
        min_amp = -min_amp; % Reprject to real angles

        max_amp = mod(max_amp, 2*pi);
        max_over_ind = max_amp > pi;
        max_amp(max_over_ind) = 2*pi - max_amp(max_over_ind);
        
        min_amp = mod(min_amp, 2*pi);
        min_over_ind = min_amp > pi;
        min_amp(min_over_ind) = 2*pi - min_amp(min_over_ind);

        
        all_amp = [max_amp ; min_amp];
        %all_amp_T = [max_amp_T ; min_amp_T];
        
        [std_amp, mean_amp] = std(all_amp);
        [std_amp_dev, mean_amp_dev] = std(abs(all_amp - params.ref_angle));

        %expFit = fit(all_amp_T, param.angle-all_amp, "exp1");

        newout.data.amp_dev = abs(all_amp - params.ref_angle);
    catch E
        % Return whats possible (avoid crashes)
        newout.val =   [mean_gsyn_pos, mean_gsyn_neg, std_gsyn_pos, std_gsyn_neg, ...
                        mean_gsm_mot_pos, mean_gsm_mot_neg, std_gsm_mot_pos, std_gsm_mot_neg, ...
                        NaN, NaN, NaN, NaN,...
                        NaN, NaN, NaN, NaN];
        return
    end

    try % Try-Catch of slip detection
        % === Find Mean and Std of Slip ===
    
        % Extract Data
        filter_index = find(o.hco_out.time >= startAnalyseTime);
        T_hco = o.hco_out.time(filter_index);
        hco = squeeze(o.hco_out.signals.values);
        hco_pos = hco(filter_index, 1);
        hco_neg = hco(filter_index, 2);
    
        filter_index = find(o.event_out.time >= startAnalyseTime);
        T_event = o.event_out.time(filter_index);
        event = squeeze(o.event_out.signals.values);
        event_pos = event(filter_index, 1);
        event_neg = event(filter_index, 2);
    
        % Find start of bursts 
        [up_crossings_hco_pos, down_crosssings_hco_pos] = ...
                compute_transition_times(T_hco,hco_pos);
        [up_crossings_hco_neg, down_crosssings_hco_neg] = ...
                compute_transition_times(T_hco,hco_neg);
    
        pos_hco_start = compute_burst_start_time(up_crossings_hco_pos, ...
                                                 down_crosssings_hco_pos);
        neg_hco_start = compute_burst_start_time(up_crossings_hco_neg, ...
                                                 down_crosssings_hco_neg);
        % The Neurons are Not Bursting => Abort Slip Detection
        if isempty(pos_hco_start) || isempty(neg_hco_start)
            error("No bursting activity in HCO neurons");
        end
    
        % Find start of synchronization events
    
        [pos_event_start, ~] = compute_transition_times(T_event,event_pos);
        [neg_event_start, ~] = compute_transition_times(T_event,event_neg);
    
        % Compute Slip

        if pos_hco_start(1) < neg_hco_start(1)
            state = 0; 
            % Remove all negative event before the (first) positive burst
            neg_event_start = neg_event_start(neg_event_start > pos_hco_start(1));
            for i = 1:length(neg_hco_start) % Check alternance
                if pos_hco_start(i) > neg_hco_start(i)
                    error("Bursts do not happen in alternance");
                end
            end
        else
            state = 1; 
            % Remove all positive event before the (first) negative burst
            pos_event_start = pos_event_start(pos_event_start > neg_hco_start(1));
            for i = 1:length(pos_hco_start) % Check alternance
                if neg_hco_start(i) > pos_hco_start(i)
                    error("Bursts do not happen in alternance");
                end
            end
        end

        if pos_hco_start(end) > neg_hco_start(end)
            % Remove all negative event after the (last) positive burst
            neg_event_start = neg_event_start(neg_event_start < pos_hco_start(end));
        else 
            % Remove all positive event after the (last) negative burst
            pos_event_start = pos_event_start(pos_event_start < neg_hco_start(end));
        end


        pos_slip = zeros(length(pos_event_start),1);
        ind_neg_hco = 0;
        ind_pos_event = 0; 
        while ind_pos_event < length(pos_event_start)
            ind_pos_event = 1+ind_pos_event;
            event_time = pos_event_start(ind_pos_event);
            % Find between which negative events we are
            while ind_neg_hco <  length(neg_hco_start) && ...
                  event_time > neg_hco_start(ind_neg_hco+1)
                ind_neg_hco =1+ind_neg_hco;
            end
            pos_slip(ind_pos_event) = pos_hco_start(ind_neg_hco+1 - state) - ...
                                   event_time;
                
        end

        neg_slip = zeros(length(neg_event_start),1);
        ind_pos_hco = 0;
        ind_neg_event = 0; 
        while ind_neg_event < length(neg_event_start)
            ind_neg_event = 1+ind_neg_event;
            event_time = neg_event_start(ind_neg_event);
            % Find between which positive events we are
            while ind_pos_hco <  length(pos_hco_start) && ...
                  event_time > pos_hco_start(ind_pos_hco+1)
                ind_pos_hco =1+ind_pos_hco;
            end
            neg_slip(ind_neg_event) = neg_hco_start(ind_pos_hco+1 + state-1) - ...
                                   event_time;
                
        end
        
        all_slip = [pos_slip ; neg_slip];
    
        [std_slip, mean_slip] = std(all_slip);
        [std_slip_abs, mean_slip_abs] = std(abs(all_slip));
    
        newout.data.freq_slip = abs(all_slip);
    catch E
        % Return whats possible (avoid crashes)
        newout.val =   [mean_gsyn_pos, mean_gsyn_neg, std_gsyn_pos, std_gsyn_neg, ...
                        mean_gsm_mot_pos, mean_gsm_mot_neg, std_gsm_mot_pos, std_gsm_mot_neg, ...
                        mean_amp, std_amp, mean_amp_dev, std_amp_dev, ...
                        NaN, NaN, NaN, NaN];
        return;
    end

    % === Export Results ===
    newout.val =   [mean_gsyn_pos, mean_gsyn_neg, std_gsyn_pos, std_gsyn_neg, ...
                    mean_gsm_mot_pos, mean_gsm_mot_neg, std_gsm_mot_pos, std_gsm_mot_neg, ...
                    mean_amp, std_amp, mean_amp_dev, std_amp_dev, ...
                    mean_slip, std_slip, mean_slip_abs, std_slip_abs];
    %newout.pos_hco_start = pos_hco_start;
    %newout.pos_event_start = pos_event_start;
    %newout.neg_hco_start = neg_hco_start;
    %newout.neg_event_start = neg_event_start;  
    %newout.pos_slip = pos_slip;
    %newout.neg_slip = neg_slip;
end


function [up_crossings, down_crossings] = compute_transition_times(T, Y_logical)
    prev = Y_logical(1);
    t_prev = T(1);

    up_crossings = zeros(length(T), 1);
    down_crossings = zeros(length(T), 1);
    ind_upcross = 0;
    ind_downcross = 0;

    for i = 2:length(T)
        curr = Y_logical(i);
        if ~prev && curr
            ind_upcross = ind_upcross + 1;
            % Linear interpolation
            %t_prev = T(i-1);
            t_curr = T(i);
            t_inter = (t_prev + t_curr) / 2;
            up_crossings(ind_upcross) = t_inter;
        elseif prev && ~curr
            ind_downcross = ind_downcross + 1;
            % Linear interpolation
            %t_prev = T(i-1);
            t_curr = T(i);
            t_inter = (t_prev + t_curr) / 2;
            down_crossings(ind_downcross) = t_inter;
        end
        t_prev = T(i);
        prev = curr;
    end
    up_crossings = up_crossings(1:ind_upcross);
    down_crossings = down_crossings(1:ind_downcross);
end


function [up_crossings, down_crossings] = compute_zero_crossings_time(T, Y)
    prev = Y(1);
    t_prev = T(1);

    up_crossings = zeros(length(T), 1);
    down_crossings = zeros(length(T), 1);
    ind_upcross = 0;
    ind_downcross = 0;

    for i = 2:length(T)
        curr = Y(i);
        if curr ~= 0 % Stall if current value is zero
            if prev < 0 && curr > 0
                ind_upcross = ind_upcross + 1;
                % Linear interpolation
                %t_prev = T(i-1);
                t_curr = T(i);
                t_inter = (t_prev*curr - t_curr*prev) / (curr - prev);
                up_crossings(ind_upcross) = t_inter;
            elseif prev > 0 && curr < 0
                ind_downcross = ind_downcross + 1;
                % Linear interpolation
                %t_prev = T(i-1);
                t_curr = T(i);
                t_inter = (t_prev*curr - t_curr*prev) / (curr - prev);
                down_crossings(ind_downcross) = t_inter;
            end
            t_prev = T(i);
            prev = curr;
        end
    end
    up_crossings = up_crossings(1:ind_upcross);
    down_crossings = down_crossings(1:ind_downcross);
end

function bursts_start = compute_burst_start_time(up_crossings, down_crossings)
    % Alternance up -> down -> up -> down ... (up -> down) or (up-> down -> up)
    if up_crossings(1) < down_crossings(1)
        down_sep = up_crossings(2:end) - down_crossings(1:length(up_crossings)-1);
    % Alternance down -> up -> down ... (down-> up -> down) or (down -> up)
    else
        down_sep = up_crossings(1:end) - down_crossings(1:length(up_crossings));
    end
    
    % Group Down time 
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
    
    % Check that we have burst
    if sum_max_down/num_max_down > 4*sum_min_down/num_min_down
        bursts_start_map = zeros(length(up_crossings),1);
        if up_crossings(1) < down_crossings(1)
        % down_sep(i) ends with up_crossings(i+1) 
        % and length(down_sep) = length(up_crossings) - 1
            bursts_start_map(idx_down(ind_max_down:end)+1) = 1;
        else
        % down_sep(i) ends with up_crossings(i)
        % and length(down_sep) = length(up_crossings)
            bursts_start_map(idx_down(ind_max_down:end)) = 1;
        end
        bursts_start_map = logical(bursts_start_map);
        bursts_start = up_crossings(bursts_start_map);
    else
        bursts_start = [];
    end
end