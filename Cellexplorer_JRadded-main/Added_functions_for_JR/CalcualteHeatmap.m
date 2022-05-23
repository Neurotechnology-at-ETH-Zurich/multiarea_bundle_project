function CalcualteHeatmap(spikes,grid_size,basename,basepath)
    
    load(fullfile(basepath,[basename,'.trace.mat']));
    
    
    x_lim = [ceil(50/grid_size),floor(550/grid_size)];
    y_lim = [ceil(100/grid_size),floor(600/grid_size)];
    
    bin_x = 0:grid_size:700;
    bin_y = 0:grid_size:700;
    hmaps = cell(1,spikes.numcells);
    hmap_size = 700/grid_size;
    

    xy_nan = find(isnan(trace_mat(:,1)));
    trace_mat(xy_nan,:) = [];
    %% generating time matrix
    time_table = zeros(hmap_size,hmap_size);
    for k = 3:length(trace_mat(:,1))
        this_x = trace_mat(k,1);
        this_y = trace_mat(k,2);
        where_x = discretize(this_x,bin_x);
        where_y = discretize(this_y,bin_y);
        if (~isnan(where_x))&&(~isnan(where_y))
        time_table(where_x,where_y) = time_table(where_x,where_y)+(trace_mat(k,3)-trace_mat(k-1,3));
        end
    end
    
    for i = 1:spikes.numcells
        xpos = spikes.x_pos{i};
        ypos = spikes.y_pos{i};
        
        blank_table = zeros(hmap_size,hmap_size);
        group_x = discretize(xpos,bin_x);
        group_y = discretize(ypos,bin_y);
        
        group_x(find(group_x>x_lim(2))) = x_lim(2);
        group_x(find(group_x<x_lim(1))) = x_lim(1);
        group_y(find(group_y>y_lim(2))) = y_lim(2);
        group_y(find(group_y<y_lim(1))) = y_lim(1);



        for j = 1:length(xpos)
            if (~isnan(group_x(j)))&&(~isnan(group_y(j)))
            blank_table(group_x(j),group_y(j)) = blank_table(group_x(j),group_y(j))+1;
            end
        end

        for k = 1:hmap_size
            for l = 1:hmap_size
                if blank_table(k,l) ~= 0
                    if time_table(k,l) == 0
                    blank_table(k,l) = 0;
                    else
                    blank_table(k,l) = blank_table(k,l)/time_table(k,l);
                    end
                end
            end
        end
           

        hmaps{i} = blank_table;
        
    end
    save(fullfile(basepath,[basename,'.hmaps.cellinfo.mat']),'hmaps');
end