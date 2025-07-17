function [ siteLoc,siteMap,ignoreSites,nChans,sampleRate,shankMap]=import_prm(file_name)

prm_file = textread(char(file_name), '%s','bufsize',40950000, 'delimiter', '\n');
%
% mPFC_channel_num = cellfun(@str2num,[sessionInfo.ElecGp{1,4}.channel]);
%

for sorok_prb=1:length(prm_file)


    shank_hely = strfind (prm_file{sorok_prb}, 'shankMap = ');
    if ~isempty (shank_hely)
        spaces_4 = findstr(prm_file{sorok_prb}, ',');
        prm_file{sorok_prb}(spaces_4)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '[');
        zarojel_stop = findstr(prm_file{sorok_prb}, ']');
        shankMap_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        shankMap  =sscanf(  shankMap_string, '%d');
    end

    shank_hely = strfind (prm_file{sorok_prb}, 'siteLoc = ');
    if ~isempty (shank_hely)
        spaces_4 = findstr(prm_file{sorok_prb}, ',');
        prm_file{sorok_prb}(spaces_4)=' ';
        spaces_1 = findstr(prm_file{sorok_prb}, ';');
        prm_file{sorok_prb}(spaces_1)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '[');
        zarojel_stop = findstr(prm_file{sorok_prb}, ']');

        siteLoc_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        siteLoc =sscanf(  siteLoc_string, '%d');
        siteLoc=(reshape(siteLoc,[2,256]))';
    end

    shank_hely = strfind (prm_file{sorok_prb}, 'siteMap = ');
    if ~isempty (shank_hely)
        spaces_4 = findstr(prm_file{sorok_prb}, ',');
        prm_file{sorok_prb}(spaces_4)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '[');
        zarojel_stop = findstr(prm_file{sorok_prb}, ']');

        siteMap_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        siteMap  =sscanf(  siteMap_string, '%d');
    end

    shank_hely = strfind (prm_file{sorok_prb}, 'nChans = ');
    if ~isempty (shank_hely)
        %         spaces_4 = findstr(prm_file{sorok_prb}, ';');
        %         prm_file{sorok_prb}(spaces_4)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '= ');
        zarojel_stop = findstr(prm_file{sorok_prb}, ';');

        nChans_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        nChans   =sscanf(nChans_string, '%d');
    end


    shank_hely = strfind (prm_file{sorok_prb}, 'sampleRate = ');
    if ~isempty (shank_hely)
        %         spaces_4 = findstr(prm_file{sorok_prb}, ';');
        %         prm_file{sorok_prb}(spaces_4)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '= ');
        zarojel_stop = findstr(prm_file{sorok_prb}, ';');

        sampleRate_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        sampleRate   =sscanf(sampleRate_string, '%d');
    end

    ignoreSites_hely = strfind (prm_file{sorok_prb}, 'ignoreSites = [');
    if ~isempty (ignoreSites_hely)
        spaces_4 = findstr(prm_file{sorok_prb}, ',');
        prm_file{sorok_prb}(spaces_4)=' ';
        zarojel_start = findstr(prm_file{sorok_prb}, '[');
        zarojel_stop = findstr(prm_file{sorok_prb}, ']');

        ignoreSites_string= prm_file{sorok_prb}(zarojel_start+1:zarojel_stop-1);
        ignoreSites =sscanf(  ignoreSites_string, '%d');
    end
end