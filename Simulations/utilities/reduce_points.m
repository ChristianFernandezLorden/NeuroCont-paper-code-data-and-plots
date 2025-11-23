% Based on Ramer-Douglas-Peucker Algorithm

function [reduced_vec, reduced_t] = reduce_points(original_vec, original_t, epsilon)
    range = max(original_vec) - min(original_vec);
    min_slope_diff = epsilon*range;
    
    reduced_values = true(1,length(original_vec));

    stack_start = java.util.Stack();
    stack_end = java.util.Stack();

    stack_start.push(1);
    stack_end.push(length(original_vec));

    while ~stack_start.empty()
        start_ind = stack_start.pop();
        end_ind = stack_end.pop();

        d_max = 0;
        index_max = 0;
        for i = start_ind+1:end_ind-1
            d = perpendicular_distance( ...
                    original_t(start_ind), original_vec(start_ind),...
                    original_t(end_ind), original_vec(end_ind),...
                    original_t(i), original_vec(i)...
                );
            if d > d_max
                d_max = d;
                index_max = i;
            end
        end

        if d_max > epsilon
            % First half
            stack_start.push(start_ind);
            stack_end.push(index_max);
            % Second half
            stack_start.push(index_max);
            stack_end.push(end_ind);
        else
            for i = start_ind+1:end_ind-1
                reduced_values(i) = false;
            end
        end
    end
    
    reduced_vec = original_vec(reduced_values);
    reduced_t = original_t(reduced_values);
end

function d = perpendicular_distance(xl1, yl1, xl2, yl2, x, y)
    d = abs((xl2-xl1)*(yl1-y) - (xl1-x)*(yl2-yl1));
    d = d/sqrt((xl2-xl1)*(xl2-xl1) + (yl2-yl1)*(yl2-yl1));
end

