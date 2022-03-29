function prmdata = readprm (filename)
% this is a function to read important data from basename.prm
% set parameter names as substrings 
  parameter_num = 4;

  shank_map = 'shankMap';
  site_loc = 'siteLoc';
  site_map = 'siteMap';
  evtWindow_line = 'evtWindow';

  prm_file = fopen(filename);
  tline = fgetl(prm_file);

  while ischar(tline)
     if contains(tline,shank_map)
        eval(tline);
     elseif contains(tline,site_loc)
        eval(tline);
     elseif contains(tline,site_map)
         eval(tline);
     elseif contains(tline,evtWindow_line)
         eval(tline);
     end
     tline = fgetl(prm_file);
  end
  
  prmdata = cell(1,parameter_num);
  prmdata{1} = shankMap;
  prmdata{2} = siteLoc;
  prmdata{3} = siteMap;
  prmdata{4} = evtWindow;
end