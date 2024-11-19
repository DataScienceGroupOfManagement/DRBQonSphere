function [Nvec, indjCell] = RandSplitDataInner(Ntr, m, N0min)

if nargin < 3  % uniform split
    Nvec = zeros(m,1);
    n = floor(Ntr/m); 
    for k = 1:m-1
        Nvec(k) = n;
    end
    Nvec(m) = Ntr-n*(m-1);

    indjCell =  cell(1,m);
    idx_rand = randperm(Ntr);
    count = 0;
    for k = 1:m
        indjCell{k} = idx_rand(count+1:count+Nvec(k));
        count = count+Nvec(k);
    end
else
    indSplit = [sort(round(rand(1,m-1)*Ntr)) Ntr];
    indjCell = cell(1,m);
    Nvec = zeros(m,1);
    indjCell{1} = 1:indSplit(1);
    Nvec(1) = length(indjCell{1});
    for jj = 1:m-1         
        indjCell{jj+1} = indSplit(jj)+1:indSplit(jj+1);
        Nvec(jj+1) = length(indjCell{jj+1});
    end
    
    indNeedAdd = find(Nvec<N0min);
    while ~isempty(indNeedAdd) 
        for kk = 1:length(indNeedAdd)
            [MaxNtr, indMaxNtr] = max(Nvec);
            indsample = indjCell{indMaxNtr};        
            Ntr_change = min(MaxNtr-N0min, N0min-Nvec(indNeedAdd(kk)));
                   
            indjCell{indNeedAdd(kk)} = [indjCell{indNeedAdd(kk)} indsample(1:Ntr_change)];
            indsample(1:Ntr_change) = [];
            Nvec(indNeedAdd(kk)) = length(indjCell{indNeedAdd(kk)});  
            indjCell{indMaxNtr} = indsample;
            Nvec(indMaxNtr) = length(indsample);
        end
        indNeedAdd = find(Nvec<N0min);
    end
end