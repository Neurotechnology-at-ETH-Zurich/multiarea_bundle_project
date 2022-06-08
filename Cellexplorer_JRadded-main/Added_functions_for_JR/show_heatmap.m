function show_heatmap(handles)
    plot3(handles.trace.trace_mat(:,1),handles.trace.trace_mat(:,2),repmat(min(min(handles.hmaps.hmaps{handles.cell_id}))-100,length(handles.trace.trace_mat(:,1)),1),'r');
    xlim([0,700])
    ylim([0,600])
    xlabel('coordinate of the physical location /mm')
    hold on
    surf(handles.x,handles.y,handles.hmaps.hmaps{handles.cell_id},'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(2)
    alpha(0.8)
    hold off
end