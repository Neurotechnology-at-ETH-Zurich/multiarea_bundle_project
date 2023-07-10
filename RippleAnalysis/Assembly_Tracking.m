
close all
for structures_NUM=1:4%length(structures_unique)
    figure(structures_NUM)

    %shift_subplot=[0 max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID))), max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*2,max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*3,max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*4]
    shift_subplot=[0 eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))']), eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))*2']), eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM))  ',1),Assembly_max_during_Ripple_ID)))*3']),eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))*4'])];
    for nume_files=1:numel(path)-1

        if ~isempty(strfind(Path_temp,'rTBY'))
            last_slesh=(strfind(Path_temp,'\'));
            saving_filename=Path_temp((strfind(Path_temp,'rTBY')):last_slesh(end)-1);
        else
            saving_filename='unknown';
        end
        s=1;

        cellId_text=eval(['(Identified_neruons(ismember(Identified_neruons(:,nume_files),[(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ')]),:))']);
        cellId_text= cellId_text';
        [~,sorthely]=sort(cellId_text(nume_files,:))
        cellId_text=cellId_text(:,sorthely);
        temp_orig_ID=eval(['[(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ')]']);

        cellId_text(size(cellId_text,1)+1,:)= eval(['temp_orig_ID(ismember([(Assembly_cellID(nume_files).'  char(structures_unique(structures_NUM)) ')],Identified_neruons(:,nume_files)))']);

        for assembly_num=1:eval(['size(AssemblyTemplates(nume_files).'  char(structures_unique(structures_NUM)) ',2)'])
            subplot_activity=subplot(length([path.sessionID]),eval(['max((arrayfun(@(s)size(s.'  char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))']),s+shift_subplot(nume_files))
            %      significant=eval(['logical(sum([((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num))), ((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))],2))']);

            %% Otsu's Methode for threshold
            significant=eval(['abs(AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num))>graythresh(abs(AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM))  '(:,assembly_num)))']);

            try
                [~,ia,~]=intersect(Identified_neruons(:,nume_files),eval(['Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) '(significant)']));
            catch
                disp('error')
            end
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '(assembly_num,:)={ia}']); %%%new
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_Original(assembly_num,:)={Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) '(significant)}']);
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_significant(assembly_num,:)=  {significant}']);
            clear ia;

            stem(find(significant),eval(['AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(significant,assembly_num)']),'r-','filled')
            ylim([-1 1])
            subtitle(num2str(assembly_num))
            hold on;
            stem(find(significant==0),eval(['AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(~significant,assembly_num)']),'b-','filled')

            ylim([-1 1])
            xlim([0 length(cellId_text)+2])
            for num_id=1:size(cellId_text,2)
                %   cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ])}

                session_id_for_label=([path.sessionID])
                if length(session_id_for_label)~=5
                    disp('This code need exactly 5 sessions, change it if you have more ERROR!!!!')
                end

                cellId_Text(num_id)={['S' num2str(session_id_for_label(1)) ' ID: ' num2str(cellId_text(1,num_id)) ' S' num2str(session_id_for_label(2)) ' ID: ' num2str(cellId_text(2,num_id)) ' S'  num2str(session_id_for_label(3))  ' ID: ' num2str(cellId_text(3,num_id)) ' S'  num2str(session_id_for_label(4))  ' ID: ' num2str(cellId_text(4,num_id)) ' S'  num2str(session_id_for_label(5))  ' ID: ' num2str(cellId_text(5,num_id))   ' Cells ID ' num2str(cellId_text(end,num_id))]};

            end
            if assembly_num==1
                set(subplot_activity,'XTick',1:eval(['length([(Cells(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ').CellID_original)])']),'XTickLabel', cellId_Text,'FontSize',16');
            else
                %             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).IL).CellID_original)])],'XTickLabel',{(Cells(Assembly_cellID(nume_files).IL).CellID_original)},'FontSize',16')
                set(subplot_activity,'XTick',1:eval(['length([(Cells(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ').CellID_original)])']),'XTickLabel',   cellId_text(end,:),'FontSize',16')
            end
            view(90,-90)
            s=s+1;
        end
        clear   cellId_text  cellId_Text
    end

    sgtitle(['Assembly Template in '  char(structures_unique(structures_NUM))  ])

    fig=gcf
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 4000 4000];
    fig.PaperSize = [4000 4000];
    saveas(fig,[char(saving_filename) '_Assembly_Template_' char(structures_unique(structures_NUM))],'svg')
    saveas(fig,[char(saving_filename) '_Assembly_Template_' char(structures_unique(structures_NUM))],'fig')
end
%
close all
%


for structures_NUM=1:length(structures_unique)

    clear intersect_neurons_assembly_ID  intersect_neurons_assembly_percent  sorted_intresct_assemblies_similarity
    for darab_assembly=1:eval(['numel(Significant_Neurons_Assemblies(1).' char(structures_unique(structures_NUM)) ')'])
        % calculate the % of common cells within an assemblymatrix
        Ref_assembly=eval(['Significant_Neurons_Assemblies(1).' char(structures_unique(structures_NUM)) '{darab_assembly}(:,:)']);
        for session_n=1:numel(Significant_Neurons_Assemblies)
            for  darab_assembly_futo=1:eval(['numel(Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) ')'])
                intersect_neurons_assembly_ID(session_n,darab_assembly_futo,darab_assembly)=eval(['{intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{darab_assembly_futo}(:,:))}'])
                intersect_neurons_assembly_percent(session_n,darab_assembly_futo,darab_assembly)=eval(['length(intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{darab_assembly_futo}(:,:)))./(length(Ref_assembly)+length(Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{darab_assembly_futo}(:,:))-length(intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{darab_assembly_futo}(:,:))))']);
            end
        end
    end
    intersect_neurons_assembly_percent(isnan(intersect_neurons_assembly_percent))=0;

    %     for i=1:size(intersect_neurons_assembly_percent,3)
    %         for session_n=1:size(intersect_neurons_assembly_percent,1)
    %             sorted_intresct_assemblies_similarity(session_n,:,i)=sort(intersect_neurons_assembly_percent(session_n,:,i),'descend')
    %         end
    %
    %     end

    sorted_intresct_assemblies_similarity=intersect_neurons_assembly_percent;

    s_valtozo=([2:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
    %cmap = colormap(lines);
    for assembly_num=1:size(sorted_intresct_assemblies_similarity,3)
        figure

        subplot((size(sorted_intresct_assemblies_similarity,1)),2,[1:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
        heatmap(sorted_intresct_assemblies_similarity(:,:,assembly_num))

        title({'Similarity of the assemblies across sessions',['refAssembly ID:' num2str(assembly_num) ' Structure: '  char(structures_unique(structures_NUM))]})

        [~,max_sim_assembly_ID]=max(intersect_neurons_assembly_percent(:,:,assembly_num)')

        for nume_files=1:size(sorted_intresct_assemblies_similarity,1)
            subplot_activity=subplot((size(sorted_intresct_assemblies_similarity,1)),2, s_valtozo(nume_files))

            %% significant stem plot if there is only one assembly/session.
            if size(sorted_intresct_assemblies_similarity,3)>1

                significant=eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_significant{max_sim_assembly_ID(nume_files),1}(:,:)]'])
                [~, assembly_order]=sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,max_sim_assembly_ID(nume_files))]']))

                if eval(['sum(AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM))  '(significant,max_sim_assembly_ID(nume_files))<0)>1'])
                    flip_szorzo=-1;
                else
                    flip_szorzo=1;
                end
                stem(find(significant(assembly_order)),sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(significant,max_sim_assembly_ID(nume_files))]'])).*flip_szorzo,'r-','filled');
                set(subplot_activity,'XTick',find(significant(assembly_order)),'XTickLabel',eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '{max_sim_assembly_ID(nume_files)}(:,:)]']),'FontSize',16');

                %               Significant_Neurons_Assemblies(nume_files).dHP_significant{max_sim_assembly_ID(nume_files),1}(:,:)
                %                 [Significant_Neurons_Assemblies(nume_files).dHP{max_sim_assembly_ID(nume_files)}(:,:)]
            else size(sorted_intresct_assemblies_similarity,3)<1 %It did not use the maximum tracakble assebmly member/session
                significant=eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_significant{size(sorted_intresct_assemblies_similarity,3),1}(:,:)]'])
                [~, assembly_order]=sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,size(sorted_intresct_assemblies_similarity,3))]']))
               
                if eval(['sum(AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM))  '(significant,size(sorted_intresct_assemblies_similarity,3))<0)>1'])
                    flip_szorzo=-1;
                else
                    flip_szorzo=1;
                end
                stem(find(significant(assembly_order)),sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(significant,size(sorted_intresct_assemblies_similarity,3))]'])).*flip_szorzo,'r-','filled');
                set(subplot_activity,'XTick',find(significant(assembly_order)),'XTickLabel',eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '{size(sorted_intresct_assemblies_similarity,3)}(:,:)]']),'FontSize',16');
            end

            %% non significant plot
            hold on;
            
            if  size(sorted_intresct_assemblies_similarity,3)>1
                [~, assembly_order]=sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,max_sim_assembly_ID(nume_files))]']))

                stem(find(significant(assembly_order)==0),(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(~significant,max_sim_assembly_ID(nume_files))]'])).*flip_szorzo,'b-','filled')
                subtitle(['Assembly ID: ' num2str(max_sim_assembly_ID(nume_files))])
            else  size(sorted_intresct_assemblies_similarity,3)<1
                [~, assembly_order]=sort(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,size(sorted_intresct_assemblies_similarity,3))]']))
               
                stem(find(significant(assembly_order)==0),(eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(~significant,size(sorted_intresct_assemblies_similarity,3))]'])).*flip_szorzo,'b-','filled')
                subtitle(['Assembly ID: ' num2str(size(sorted_intresct_assemblies_similarity,3))])
            end

            xlim([0 length(significant)+1])
            ylim([-1 +1])

            fig_ax=gca;
            fig_ax.PlotBoxAspectRatio=([3,2,3]) ;

            view(90,-90)
        end
        fig=gcf
        fig.PaperUnits = 'points';
        fig.PaperPosition = [0 0 1400 1000];
        fig.PaperSize = [1400 1000];
        saveas(fig,['TBY37_Assembly_Tracking_'  char(structures_unique(structures_NUM)) '_assembly_' num2str(assembly_num)],'svg')
        saveas(fig,['TBY37_Assembly_Tracking_'  char(structures_unique(structures_NUM)) '_assembly_' num2str(assembly_num)],'tif')

    end
end

