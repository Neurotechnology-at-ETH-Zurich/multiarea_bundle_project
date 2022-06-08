function get_shank_region(handles)
    shankID = handles.shank_id;
    shank = shankID(handles.cell_id);
    switch shank
        case 1
            set(handles.shankregion,'string','iHP intermediate hippocampus');
        case 2
            set(handles.shankregion,'string','dHP dorsal hippocampus');
        case 3
            set(handles.shankregion,'string','RSC retrosplenial cortex');
        case 4
            set(handles.shankregion,'string','mPFC medial prefrontal cortex');
    end

end