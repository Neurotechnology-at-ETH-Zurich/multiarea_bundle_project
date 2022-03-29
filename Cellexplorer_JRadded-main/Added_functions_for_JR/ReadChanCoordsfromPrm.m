function ReadChanCoordsfromPrm (session)
%this is a function called in sessionTemplate.m where channel coordinates
%are generated from local prm file
prmdata = readprm([session.general.name,'.prm']);
siteLoc = prmdata{2};

chanCoords = struct;
chanCoords.x = siteLoc(:,1);
chanCoords.y = siteLoc(:,2);
chanCoords.source = [session.general.name,'.prm'];
chanCoords.layout = 'staggered';
chanCoords.shankSpacing = 1000;
delete(fullfile(session.general.basePath,[session.general.name,'.chanCoords.channelInfo.mat']));
save(fullfile(session.general.basePath,[session.general.name,'.chanCoords.channelInfo.mat']),'chanCoords');
end