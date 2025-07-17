(sum([Cells(find(ismember({Cells.Structre},'IL'))).clusterID_kmeans]==1)./sum([ismember({Cells.Structre},'IL')])).*100   


(sum([Cells(find(ismember({Cells.Structre},'IL'))).clusterID_kmeans]==1)./sum([ismember({Cells.Structre},'IL')])).*100   

lengths = sum(unique(arrayfun(@(x) size(Cells(x).individual_trials,1), 1:numel(Cells))))


(sum([Cells(find(ismember({Cells.Structre},'PrL'))).clusterID_kmeans]==1)./sum([ismember({Cells.Structre},'PrL')])).*100   


T = struct2table(statistic.RSC_dHP.mean,)



(sum([Cells(find(ismember({Cells.Structre},'IL'))).clusterID_kmeans]==1)./sum([ismember({Cells.Structre},'IL')])).*100   