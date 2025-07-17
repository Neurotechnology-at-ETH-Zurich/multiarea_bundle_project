%% Check platform

if ismac
    slesh=('/');
elseif isunix
    slesh=('/');
elseif ispc
    slesh=('\');
else
    disp('Platform not supported')
end

%% Loading variabals

[file,base_folder,indx] = uigetfile('*.mat');
listing=dir(base_folder);
cd(base_folder);

disp('Loading files............')
% load(file, 'Activities')
% load(file, 'Activitymatrix')
% load(file, 'AssemblyTemplates')
% load(file, 'Assembly_activity')
% load(file, 'Assembly_cellID')
% load(file, 'Assembly_decimal_ID')
% load(file, 'Assembly_decimal_rand_ID')
% load(file, 'Assembly_max_during_Ripple_ID')
% load(file, 'Assembly_max_random_ID')
% load(file, 'Identified_Neurons')
% load(file, 'Identified_neruons')
% load(file, 'edges_CCA')
% load(file, 'Cells')
% load(file, 'Ripples')

% getting the folder with the data, eliminate noise from path
path_counter=1;
for i=1:numel(listing)
    k = strfind(listing(i).name,'freely_behav');
    if ~isempty(k)
        path(path_counter).folders=fullfile(base_folder,listing(i).name);
        l=strfind(listing(i).name,'_');
        path(path_counter).sessionID=str2num(listing(i).name(1:l(1)-1));
        path_counter=path_counter+1;
    end

end
clear path_counter;

% ordered of the session by session ID
[~,sorting_session]=sort([path.sessionID]);
path=path(sorting_session);

%% make a folder for the figures

dirFlags = [listing.isdir];
sum(ismember({listing(dirFlags).name},'Figures'));
if sum(ismember({listing(dirFlags).name},'Figures'))~=1
    [status, msg, msgID] = mkdir(['Figures']); %  slesh 'Assembly_activity'
end

%% Constant:
edges_for_assembly=edges_CCA; %detection windows around the SPW-Rs Inherited from Assembly_merged.m
start_stop_rip=[find(round(edges_for_assembly,2)==-0.13) find(round(edges_for_assembly,2)==0.13)];
STAT={};
flag_plot=1; %plot or not
len=length(edges_for_assembly)*2;

structures_unique=unique({Cells.Structre}); %get all name for structure within Cells
color_lines=colormap('lines');

%% Start Analysis
for nume_files=1:numel(path)

    Path_temp=char(path(nume_files).folders);

    %platform /\ check and create filename -> for savings
    if ~isempty(strfind(Path_temp,'rTBY'))
        if ismac
            last_slesh=(strfind(Path_temp,'/'));
        elseif isunix
            last_slesh=(strfind(Path_temp,'/'));
        elseif ispc
            last_slesh=(strfind(Path_temp,'\'));
        else
            disp('Platform not supported');
        end
        saving_filename=Path_temp((strfind(Path_temp,'rTBY')):last_slesh(end)-1);
    end

    for structures_NUM=1:length(structures_unique)
        disp(['Working on Session ID: ' num2str(path(nume_files).sessionID) ' Sturcture name: ' char(structures_unique(structures_NUM)) ' - Threshold is 95th of baseline']);
        P=[];
        if ~isempty(structures_unique(structures_NUM))

            for assembly_num_total=1:eval(['(size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3))'])
      trashold=eval(['2*mean(((std(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM})  '((len/2):length((start_stop_rip(1):start_stop_rip(2)))+(len/2),:,assembly_num_total))))));']);
  
      %  trashold=eval(['prctile(mean(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(:,:,assembly_num_total)),95);']);
      %  trashold=eval(['prctile(mean(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '((len/2):length((start_stop_rip(1):start_stop_rip(2)))+(len/2),:,assembly_num_total)),95);']);
       
                eventSubset_temp(assembly_num_total,:)=eval(['(mean(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '((start_stop_rip(1):start_stop_rip(2))+(len/2),:,assembly_num_total))))>trashold;']);
                eventSubset_temp_before_binary(assembly_num_total,:)=eval(['(mean(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(start_stop_rip(1):start_stop_rip(2)),:,assembly_num_total))))>trashold;']);
                eventSubset_max(assembly_num_total,:)=eval(['(mean(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '((start_stop_rip(1):start_stop_rip(2))+(len/2),:,assembly_num_total))));']);
                eventSubset_max_before(assembly_num_total,:)=eval(['(mean(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(start_stop_rip(1):start_stop_rip(2)),:,assembly_num_total))));']);

            end
            minDigits=size(eventSubset_temp,1);

            if eval(['(size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3))'])>1
                temp_binary=sum(repmat((pow2(size(eventSubset_temp,1)-1:-1:0))',[1,size(eventSubset_temp,2)]).*eventSubset_temp);
                temp_binary_before=sum(repmat((pow2(size(eventSubset_temp_before_binary,1)-1:-1:0))',[1,size(eventSubset_temp_before_binary,2)]).*eventSubset_temp_before_binary);
            else
                temp_binary=(repmat((pow2(size(eventSubset_temp,1)-1:-1:0))',[1,size(eventSubset_temp,2)]).*eventSubset_temp);
                temp_binary_before=(repmat((pow2(size(eventSubset_temp_before_binary,1)-1:-1:0))',[1,size(eventSubset_temp_before_binary,2)]).*eventSubset_temp_before_binary);
            end

            eval(['Assembly_decimal_ID(nume_files).' char(structures_unique{structures_NUM}) '= temp_binary;']);
            eval(['Assembly_decimal_rand_ID(nume_files).' char(structures_unique{structures_NUM}) '= temp_binary_before;']);

            eval(['Assembly_max_during_Ripple_ID(nume_files).' char(structures_unique{structures_NUM}) '= eventSubset_max;']);
            eval(['Assembly_max_random_ID(nume_files).' char(structures_unique{structures_NUM}) '=  eventSubset_max_before;']);

            if sum(temp_binary)>0
                trash=10;
                [N_peri,edges] = histcounts(temp_binary,0.5:1:max(temp_binary)+0.5);
                edges_peri=floor(edges(2:end));
                [N_pre,edges] = histcounts(temp_binary_before,0.5:1:max(temp_binary)+0.5);
                edges_pre=floor(edges(2:end));
                binary_num_peri=find(and(N_peri>trash,N_peri>N_pre));
                extra_assembly_during_ripple=(setxor(find(N_pre),find(N_peri)));
                extra_assembly_during_ripple=extra_assembly_during_ripple(N_peri(setxor(find(N_pre),find(N_peri)))>0);
                extra_assembly_during_ripple=extra_assembly_during_ripple(N_peri(extra_assembly_during_ripple)>trash);

                binary_num_peri=[binary_num_peri extra_assembly_during_ripple];

                %% PLotting the number of activated assembly events by specific Ripples
                if flag_plot==1
                    fig_combination=figure;
                    elotte= histcounts(temp_binary_before(temp_binary_before>0),0.5:1:max(temp_binary)+0.5);
                    hold on;
                    alatta= histcounts(temp_binary(temp_binary>0),0.5:1:max(temp_binary)+0.5);
                    hold on;
                    try
                        bar(1:length([elotte; alatta])/2, [elotte; alatta]');
                    catch
                        bar([elotte; alatta]');
                    end

                    try
                        xticks(unique(sort(binary_num_peri)))
                    catch
                        disp('error')
                    end
                    xticklabels({dec2bin(unique(sort(binary_num_peri)),minDigits)} );
                    xtickangle(90);
                    ylabel('Number of ripples','FontSize',16);
                    xlabel('Combination of selectively activited assemblies patterns before/during SWRs','FontSize',16);
                    title(char(structures_unique{structures_NUM}),'FontSize',18)
                    set(gca,'FontSize',16)
                    axis square
                    fig=gcf;

                    fig.PaperUnits = 'points';
                    fig.Renderer='painters'
                    fig.PaperPosition = [0 0 800 800];
                    fig.PaperSize = [800 800];
                    saveas(fig,[char(saving_filename) '_' num2str(path(nume_files).sessionID) '_Assembly_Probabilty_' char(structures_unique{structures_NUM})],'svg');
                    saveas(fig,[ char(saving_filename) '_' num2str(path(nume_files).sessionID) '_Assembly_Probabilty_' char(structures_unique{structures_NUM})],'tif');
                end


                %% Plot Assembly Strength

                text_for_binary=dec2bin(binary_num_peri,minDigits);
                for i= 1:length(binary_num_peri)
                    clear  Total_Ripple_selected_assembly   Total_Ripple_selected_assembly_random   Ripples_selected_for_assembly_random Ripples_selected_for_assembly

                    ripple_selected=find(ismember((Ripples(nume_files).ripples.eventID),{'bzHigh_dHP','bzLow_dHP','bzHigh_iHP','bzLow_iHP'}));
                    eventSubset=intersect(ripple_selected,(find(temp_binary==binary_num_peri(i))));
                    eventSubset_rand_temp=intersect(ripple_selected,(find(temp_binary~=binary_num_peri(i))));

                    % sum(~ismember(eventSubset,(Ripples(nume_files).ripples.common_ripples(:,1))))
                    % Ripples(nume_files).ripples.eventID((intersect(eventSubset,(Ripples(nume_files).ripples.common_ripples(:,1)))))
                    % (intersect(eventSubset,(Ripples(nume_files).ripples.common_ripples(:,2))))

                    if flag_plot==1
                        figure(i)
                        subp_ax1=subplot(2,3,1);

                        try
                            bar(1:length([elotte; alatta])/2, [elotte; alatta]');
                            ylim([0 round(max([elotte, alatta]),2,"significant")]);
                        catch
                            bar([elotte; alatta]');
                            ylim([0 round(max([elotte, alatta]),2,"significant")]);
                        end


                        try
                            xticks(unique(sort(binary_num_peri)));
                        catch
                            disp('error')
                        end
                        xticklabels({dec2bin(unique(sort(binary_num_peri)),minDigits)} );
                        xtickangle(90);
                        ylabel('Number of ripples','FontSize',16);
                        xlabel('Combination of selectively activited assemblies patterns before/during SWRs','FontSize',16);
                        title(char(structures_unique{structures_NUM}),'FontSize',18);
                        set(gca,'FontSize',16);
                        axis square;

                        for num_as=1:eval(['size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3)'])
                            temp_aerobar_random=eval(['Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(edges_for_assembly),eventSubset,num_as);']);
                            temp_aerobar_random_save(:,:,num_as)= temp_aerobar_random;
                            temp_aerobar_random=[];
                            axis square
                        end

                        % xlim([find(round(edges_for_assembly.*1000)==-500) find(round(edges_for_assembly.*1000)==500)])
                        %
                        % xticks([find(round(edges_for_assembly.*1000)==-500), find(round(edges_for_assembly.*1000)==-250), round(length(edges_for_assembly)/2),find(round(edges_for_assembly.*1000)==250), find(round(edges_for_assembly.*1000)==500)]);
                        % xticklabels([-500 -250,0,250,500])
                        % xlabel('Time [ms]')
                        % ylabel('Avaraged Assembly expression strength')

                        subp_ax2=subplot(2,3,2);
                        for num_as=1:eval(['size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3)'])
                            %hold on; plot(edges_for_assembly,mean(assembly_activity_IL(:,eventSubset_IL,num_as),2),'-','LineWidth',2)
                            options.color_area = color_lines(num_as,:);% [243 169 114]./255;    % Orange theme
                            options.color_line = color_lines(num_as,:);%[236 112  22]./255;
                            options.error='std';
                            options.alpha      = 0.5;
                            options.line_width = 2;
                            options.error      = 'std';
                            temp_aerobar=eval(['Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(length(edges_for_assembly)+1:end,eventSubset,num_as);']);
                            plot_areaerrorbar(temp_aerobar',options);

                            hold on; plot(mean(squeeze(mean(temp_aerobar_random_save,2)),2),'--k','LineWidth',2)
                            std_temp=2*mean(std(squeeze(mean(temp_aerobar_random_save,2))));
                            hold on; plot(mean(squeeze(mean(temp_aerobar_random_save,2)),2)+ std_temp,'--','LineWidth',2,'Color',  [0.800000011920929 0.800000011920929 0.800000011920929]);
                            hold on; plot(mean(squeeze(mean(temp_aerobar_random_save,2)),2)- std_temp,'--','LineWidth',2,'Color',  [0.800000011920929 0.800000011920929 0.800000011920929]);

                            temp_aerobar=[];

                            % ylim([subp.YLim(1) subp.YLim(2)]);
                            ylim([-0.5 11])
                            axis square;
                        end
                        clear   temp_aerobar_random_save
                        %linkaxes([subp_ax2, subp_ax1],'y');

                        xlim([find(round(edges_for_assembly.*1000)==-500) find(round(edges_for_assembly.*1000)==500)]);
                        xticks([find(round(edges_for_assembly.*1000)==-500), find(round(edges_for_assembly.*1000)==-250), round(length(edges_for_assembly)/2),find(round(edges_for_assembly.*1000)==250), find(round(edges_for_assembly.*1000)==500)]);
                        xticklabels([-500 -250,0,250,500]);
                        xlabel('Time (ms)');
                        %           ylabel('Assembly expression strength')
                        %
                        subplot(2,3,3)

                        %subplot(1,2,1)
                        %mean(eventSubset_max)

                        % figure;

                        DATA_for_Violion=[];
                        VariableNames={};
                        groupcol={};
                        groupcol_struc={};
                        for structures_NUM_temp=1:length(structures_unique)
                            eventSubset_max_test=[];
                            eventSubset_max_before_test=[];
                            for assembly_num_total=1:eval(['(size(Assembly_activity(nume_files). ' char(structures_unique{structures_NUM_temp}) ',3))'])

                                eventSubset_max_test(assembly_num_total,:)=eval(['(mean(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM_temp}) '((start_stop_rip(1):start_stop_rip(2))+(len/2),:,assembly_num_total))));']);
                                VariableNames=[VariableNames,[structures_unique{structures_NUM_temp} ' Assembly ' num2str( assembly_num_total)]];
                                %     VariableNames=[VariableNames,[structures_unique{structures_NUM_temp} ' BaseLine Assembly ' num2str( assembly_num_total)]];
                                DATA_for_Violion=[DATA_for_Violion,(eventSubset_max_test(assembly_num_total,temp_binary==binary_num_peri(i)))'];
                                %    DATA_for_Violion=[DATA_for_Violion,(eventSubset_max_before(assembly_num_total,temp_binary==binary_num_peri(i)))'];
                                groupcol =  [groupcol, repmat({[structures_unique{structures_NUM_temp} ' Assembly ' num2str( assembly_num_total)]},1,length((eventSubset_max_test(assembly_num_total,temp_binary==binary_num_peri(i)))'))];
                                %   groupcol =  [groupcol, repmat({[structures_unique{structures_NUM_temp} ' BaseLine Assembly ' num2str( assembly_num_total)]},1,length((eventSubset_max_test(assembly_num_total,temp_binary==binary_num_peri(i)))'))]
                                groupcol_struc=[groupcol_struc,{structures_unique{structures_NUM_temp}}];
                                %  groupcol_struc=[groupcol_struc,{structures_unique{structures_NUM_temp}}];

                                if ismember(structures_unique(structures_NUM_temp),structures_unique(structures_NUM))
                                    VariableNames=[VariableNames,[structures_unique{structures_NUM_temp} ' BaseLine Assembly ' num2str( assembly_num_total)]];
                                    DATA_for_Violion=[DATA_for_Violion,(eventSubset_max_before(assembly_num_total,temp_binary==binary_num_peri(i)))'];
                                    groupcol =  [groupcol, repmat({[structures_unique{structures_NUM_temp} ' BaseLine Assembly ' num2str( assembly_num_total)]},1,length((eventSubset_max_test(assembly_num_total,temp_binary==binary_num_peri(i)))'))]
                                    groupcol_struc=[groupcol_struc,{structures_unique{structures_NUM_temp}}];
                                end



                            end
                        end

                        % VariableNames=[VariableNames,repmat(['BaseLine ' char(structures_unique(structures_NUM))],[1 size(eventSubset_max_before(:,(temp_binary==binary_num_peri(i)))',2)])];
                        % DATA_for_Violion=[DATA_for_Violion,eventSubset_max_before(:,(temp_binary==binary_num_peri(i)))'];
                        % groupcol =  [groupcol, repmat({['BaseLine ' char(structures_unique(structures_NUM))]},1,length((eventSubset_max_test(assembly_num_total,temp_binary==binary_num_peri(i)))'))];
                        % groupcol_struc=[groupcol_struc,{['BaseLine ' char(structures_unique(structures_NUM))]}];


                        a2t = array2table(DATA_for_Violion,"VariableNames",  VariableNames);
                        groupcol_cat=categorical(groupcol);




                        [~,~,category_for_color]=unique(groupcol_struc);

                        color_violin=colormap('lines');
                        color_violin=color_violin(category_for_color,:);

                        vectorised_grouped_data=reshape(DATA_for_Violion,1,[])';

                        b=boxchart(groupcol_cat,vectorised_grouped_data,'GroupByColor', groupcol_cat,'BoxWidth',10,'MarkerStyle','none');
                        ylim([0 prctile(vectorised_grouped_data,99)+1]);

                         [p,~,~] = anova1(vectorised_grouped_data,groupcol_cat,'off');
                        % [c,~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni");

                        % [p,table,stats] = kruskalwallis(vectorised_grouped_data,groupcol_cat,'off')
                        % % [c,~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni");
                        % P(i,:)=p;
                        % STAT(i,:)={stats};
                        % TABLE(i,:)={table};


                        for colors=1:numel(b)
                            b(colors).BoxFaceColor=color_violin(colors,:);
                            b(colors).BoxEdgeColor=color_violin(colors,:);
                            b(colors).MarkerColor=color_violin(colors,:);
                        end
                        ylim([0 7]);
                        axis square
                    title(['Assmebies pattern: ' (text_for_binary(i,:)) ' in ' char(structures_unique{structures_NUM}) 'ANOVA p-value: ' num2str(round(p,4))])
                    ylabel('max abs Assembly expression strength')
                    xlabel('Assembly ID')
                   
               
                        match=[];
                        matchC={};
%str2num( text_for_binary)
% [~,helze_max]=max(alatta)
% ismember(binary_num_peri,helze_max)
                        matchC =strfind(a2t.Properties.VariableNames, structures_unique{structures_NUM} );
                        match  = ~cellfun('isempty', matchC);

                        subplot(2,3,4)
                        violinplot(a2t(:,[a2t.Properties.VariableNames(match)]),[a2t.Properties.VariableNames(match)],'ShowMean', true, 'ShowBox',false,'ShowMedian',false);

                        % x=vectorised_grouped_data(ismember(groupcol_cat,[a2t.Properties.VariableNames(match)]))
                        % y=groupcol_cat(ismember(groupcol_cat,[a2t.Properties.VariableNames(match)]))'

                        %        [p,anovatab,stats]= anova1(x,y)

                        [p,anovatab,stats]= anova1(vectorised_grouped_data(ismember(groupcol_cat,[a2t.Properties.VariableNames(match)])),groupcol_cat(ismember(groupcol_cat,[a2t.Properties.VariableNames(match)]))','off')

                        % [c,~,~,gnames] = multcompare(stats);
                        P(i,:)=p;
                        STAT(i,:)={stats};
                        TABLE(i,:)={table};


                        % violinplot([(eventSubset_max(:,(temp_binary==binary_num_peri(i))))', eventSubset_max_before(:,(temp_binary==binary_num_peri(i)))']);
                        %

                        title(['Assmebies pattern: ' (text_for_binary(i,:)) ' in ' char(structures_unique{structures_NUM}) ' p-value: ' num2str(round(p,4))]);
                        ylabel('max abs Assembly expression strength');
                        xlabel('Assembly ID');
                        % Value of the assembly strength ripple
                        %  subplot(1,2,2)
                        %  violinplot((eventSubset_max(:,find(temp_binary==binary_num(i-2))))')
                        % title(['Number of ripples: ' num2str(numRipple(i-2))])
                        % xlabel('Assembly ID')
                        %  subplot(2,2,3)
                        axis square

                        % subplot(2,3,4)
                        % yyaxis left
                        % bar(scores(idx_features))
                        % hold on;
                        % xlabel('Predictor rank')
                        % ylabel('Predictor importance score')
                        % xticks([1:length(idx_features)])
                        % xlim([0.5 10.5])
                        % xticklabels(strrep(Total_Ripple_selected_assembly.Properties.VariableNames(idx_features),'_','\_'))
                        % xtickangle(45)
                        %
                        % yyaxis right
                        % plot([1:length(idx_features)],(exp(-scores(idx_features))),'LineWidth',4)
                        % hold on; xline(0.05)
                        % ylabel('p-value')
                        % subtitle('Univariate feature ranking for classification using chi-square tests')
                        %                 if ~isempty((find(exp(-scores(idx))<0.05)))
                        %                     hline(0.05)
                        %                 line(find(exp(-scores(idx))<0.05))
                        %axis square

                        subplot(2,3,6)


                        %%correlation
                        % principle_corr_column=structures_NUM+num_as-i;
                        % DATA_for_Violion_corr=(DATA_for_Violion(:,~ismember(1:size(DATA_for_Violion,2),principle_corr_column)));
                        % VariableNames(~ismember(1:size(DATA_for_Violion,2),principle_corr_column))
                        %
                        % for num_assembly_total=1:size(DATA_for_Violion,2)-1
                        %     principle_corr_vector=(DATA_for_Violion(:,ismember(1:size(DATA_for_Violion,2),principle_corr_column)));
                        %
                        %     [R,P] = corrcoef(principle_corr_vector,DATA_for_Violion_corr(:, num_assembly_total));
                        %     R_value(num_assembly_total,:)=R(2)
                        %     if abs(R(2))>0.8
                        %         hold on;  scatter(principle_corr_vector,DATA_for_Violion_corr(:, num_assembly_total),'.r')
                        %         b1 = principle_corr_vector\DATA_for_Violion_corr(:, num_assembly_total)
                        %         yCalc1 = b1*principle_corr_vector;
                        %         hold on;   plot(principle_corr_vector, yCalc1,'Color','r');
                        %         title(['R:' num2str(R(2))])
                        %     else
                        %         hold on;  scatter(principle_corr_vector,DATA_for_Violion_corr(:, num_assembly_total),'.k')
                        %         b1 = principle_corr_vector\DATA_for_Violion_corr(:, num_assembly_total)
                        %         yCalc1 = b1*principle_corr_vector;
                        %         hold on;   plot(principle_corr_vector, yCalc1,'Color','k');
                        %     end
                        % end
                        % axis square
                        % ylim([0 max(max(DATA_for_Violion_corr))])
                        % xlim([min(principle_corr_vector)-0.1 max(principle_corr_vector)+0.1])

                        name_by_structure=Ripples(nume_files).ripples.eventID(eventSubset);

                        % intersect(Ripples(nume_files).ripples.common_ripples(:,:),eventSubset)
                        % Ripples(nume_files).ripples.eventID(Ripples(nume_files).ripples.common_ripples(:,1))

                        %               sum(ismember(Ripples(nume_files).ripples.common_ripples(:,2),eventSubset))+sum(ismember(Ripples(nume_files).ripples.common_ripples(:,1),eventSubset))
                        % name_by_structure(ismember(name_by_structure,{'bzLow_iHP','bzHigh_iHP'}))='iHP';
                        % name_by_structure(ismember(name_by_structure,{'bzLow_dHP','bzHigh_dHP'}))='dHP';
                        %
                        % name_by_structure;
                        % name_by_structure;


                        name_by_structure_cell=cellstr(name_by_structure);
                        name_by_structure_cell=regexprep(name_by_structure_cell,'_', ' ');
                        name_by_structure_cell=regexprep(name_by_structure_cell,'bz', '');


                        pie(categorical(regexprep(name_by_structure_cell,'_', ' ')));
                        subtitle(['Percentage of Ripples n:' num2str(length(name_by_structure)) '/' num2str(length(Ripples(nume_files).ripples.eventID))])
                        name_by_structure_cell=[];
                        % if i==length(binary_num_peri)
                        %     subplot(2,3,5)
                        %     gru=[];fold=[];
                        %     [gru,fold]=size(validationAccuracy_fold);
                        %     Group_Label=[];
                        %     for darab_gru=1:gru
                        %         Group_Label=[Group_Label  repmat({text_for_binary(darab_gru,:)},1,fold)];
                        %     end
                        %     Group_Label=categorical(Group_Label);
                        %     data_boxplot=[];
                        %     data_boxplot=reshape(validationAccuracy_fold',[],1);
                        %     %                 boxchart(validationAccuracy_fold')
                        %     boxchart(Group_Label,data_boxplot)
                        %     title ('Accuracy');
                        % end

                        %                     subplot(2,3,5)
                        %                     gru=[];fold=[];
                        %                     [gru,fold]=size(validationAccuracy_fold(i,:));
                        %                     Group_Label=[];
                        %                     for darab_gru=1:gru
                        %                         Group_Label=[Group_Label  repmat({text_for_binary(darab_gru,:)},1,fold)];
                        %                     end
                        %                     gru=[];fold=[];
                        %                     [gru,fold]=size(validationAccuracy_fold_rand(i,:));
                        %                     Group_Label_rand=[];
                        %                     for darab_gru=1:gru
                        %                         Group_Label_rand=[Group_Label_rand  repmat({['shuffled -' text_for_binary(darab_gru,:)]},1,fold)];
                        %                     end
                        %
                        %                     Group_Label_rand=categorical(Group_Label_rand);
                        %                     data_boxplot=[];
                        %                     data_boxplot=reshape(validationAccuracy_fold(i,:)',[],1);
                        %                     data_boxplot_rand=[];
                        %                     data_boxplot_rand=reshape(validationAccuracy_fold_rand(i,:)',[],1);
                        %                     %                 boxchart(validationAccuracy_fold')
                        %                     boxchart([Group_Label Group_Label_rand]' ,[data_boxplot; data_boxplot_rand])
                        %                     title (['Accuracy' ]);
                        %
                        %                     [h,p,ci,stats] = ttest(validationAccuracy_fold(i,:), validationAccuracy_fold_rand(i,:));
                        %                        title (['Accuracy p-value = ' num2str(round(p,3))]);
                        %                        ylim([0 1])
                        % validationAccuracy_fold_Total.valid(structures_NUM, nume_files)={validationAccuracy_fold(i,:)};
                        % % validationAccuracy_fold_Total.random(structures_NUM, nume_files)={validationAccuracy_fold_rand(i,:)};
                        % %                     validationAccuracy_fold_Total.ttest_p(structures_NUM, nume_files)=p;
                        %                     validationAccuracy_fold_Total.ttest_h(structures_NUM, nume_files)=h;
                        %                     validationAccuracy_fold_Total.ttest_ci(structures_NUM, nume_files)={ci};
                        %                     validationAccuracy_fold_Total.ttest_stats(structures_NUM, nume_files)={stats};

                        %         catch
                        %             disp('no data for the Accuracy')
                        % end

                        fig=gcf;
                        fig.PaperUnits = 'points';
                        fig.Renderer='painters';
                        fig.PaperPosition = [0 0 800 600];
                        fig.PaperSize = [800 600];
                        base_folder
                        saveas(fig,fullfile(base_folder,'Figures',[char(saving_filename) '_' num2str(path(nume_files).sessionID) '_'  char(structures_unique{structures_NUM}) '_Assembly_pattern_new_' num2str(dec2bin(binary_num_peri(i)))]),'fig')
                        saveas(fig,fullfile(base_folder,'Figures',[char(saving_filename) '_' num2str(path(nume_files).sessionID) '_'  char(structures_unique{structures_NUM}) '_Assembly_pattern_new_' num2str(dec2bin(binary_num_peri(i)))]),'svg')
                        saveas(fig,fullfile(base_folder,'Figures',[char(saving_filename) '_' num2str(path(nume_files).sessionID) '_'  char(structures_unique{structures_NUM}) '_Assembly_pattern_new_' text_for_binary(i,:)]),'tif')

                    end

                    clear Group_Label  a2t  Group_Label_rand eventSubset_rand_temp Ripples_selected_for_assembly entry eventSubset_rand Ripples_selected_for_assembly_random Total_Ripple_selected_assembly  Response
                end
                % STAT_TOTAL(structures_NUM, nume_files)={STAT};

                if ~isempty(P) & ~isempty(STAT)
                    eval(['statistic.' [char(saving_filename) '_' char(structures_unique(structures_NUM)) '_' num2str(path(nume_files).sessionID)] '.p=P;'])

                    % eval(['statistic.' [char(structures_unique(structures_NUM)) '_' char(saving_filename)] '.ci(num_cluster,1:2)=ci'])
                    eval(['statistic.' [char(saving_filename) '_' char(structures_unique(structures_NUM)) '_' num2str(path(nume_files).sessionID)] '.stats={STAT};'])
                    eval(['statistic.' [char(saving_filename) '_' char(structures_unique(structures_NUM)) '_' num2str(path(nume_files).sessionID)] '.table={TABLE};'])

                end
            end
        end
        close all
        clear    Response Response2 validationAccuracy_fold eventSubset_temp  eventSubset_temp_before_binary   eventSubset_max  eventSubset_max_before p P STAT stats

    end

end
