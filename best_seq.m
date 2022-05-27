function [sequence,pairs] = best_seq(sequence,pairs,recFactor,maxIm)
if isempty(pairs)% when pairs end
elseif length(pairs)<recFactor
elseif isempty(sequence)% initializing sequence
    pairs=pairs(randperm(length(pairs)),:);
    sequence(1,:)=pairs(1,:);
    pairs(1,:)=[];
    [sequence,pairs] = best_seq(sequence,pairs,recFactor,maxIm);
elseif size(pairs,1)<1
    sequence=[sequence; pairs(1,:)];
else
    ind=best_pair(sequence,pairs,maxIm);
    sequence=[sequence; pairs(ind,:)];
    pairs(ind,:)=[];
    [sequence,pairs] = best_seq(sequence,pairs,recFactor,maxIm);
end
end




function ind=best_pair(sequence,pairs,maxIm)
allseq=sequence(:);
uniq_pairs=unique(pairs(:));
allseq=[allseq; maxIm];
allseq=tabulate(allseq);
allseq(maxIm,2)=allseq(maxIm,2)-1;
onemat=allseq(find(ismember(allseq(:,1),uniq_pairs)),:);
freq=min(onemat(:,2));
pos=find(onemat(:,2)==freq);
pos=pos(randperm(length(pos)));
pos=pos(1);
value1=onemat(pos,1);
newmat=[pairs(find(pairs(:,1)==value1),:);pairs(find(pairs(:,2)==value1),:)];
uniq_pairs=unique(newmat(:));
uniq_pairs=uniq_pairs(uniq_pairs~=value1);
onemat=allseq(find(ismember(allseq(:,1),uniq_pairs)),:);
freq=min(onemat(:,2));
pos=find(onemat(:,2)==freq);
pos=pos(randperm(length(pos)));
pos=pos(1);
value2=onemat(pos,1);


if sum(ismember(pairs,[value1 value2],'rows'))
    ind=find(ismember(pairs,[value1 value2],'rows'));
elseif sum(ismember(pairs,[value2 value1],'rows'))
    ind=find(ismember(pairs,[value2 value1],'rows'));
else
    error('not a normal pair')
end
end



