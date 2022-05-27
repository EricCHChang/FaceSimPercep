function tpOrder(subID,studyID) 
% Generate order mat file for each subject in the perception-based 
% similarity rating task
% Inputs:
%   subID - a subject's ID (a string; e.g. 'ips01')
%   study - the study ID (a string; e.g., 'ips')

%%
rng('shuffle') 
% rand('state', sum(100*clock)); % Initialize the random number generator

scriptName = mfilename;
p = mfilename('fullpath');
root_dir = p(1:end-length(scriptName));
output_dir = fullfile(root_dir,['order_tp_' studyID]);

crit = 1;
while crit
    [sequence,pairs] = best_seq([],nchoosek(1:57,2),1500,57);
    [sequence,pairs] = best_seq(sequence,pairs,1100,57);
    [sequence,pairs] = best_seq(sequence,pairs,700,57);
    [sequence,pairs] = best_seq(sequence,pairs,300,57);
    [sequence,pairs] = best_seq(sequence,pairs,0,57);
    length(pairs)
    if length(pairs)==0
        crit=0;
    end
end

stimNum = 57;
blockNum = 14;
blocklength = 114;

% Make the images present roughly equal times on the left and right sides
pair = sequence;
%pair_ori = pair;
for i = 1:blockNum
    blkind = (i-1)*blocklength+1:i*blocklength;
    pair_old = pair(blkind,:);
    start = 1;
    tic
    while start ~=0
        for j = 1:stimNum
            leftNum = nnz(find(pair(blkind,1)==j));
            rightNum = nnz(find(pair(blkind,2)==j));
            if abs(leftNum-rightNum)>1
                if leftNum > rightNum
                    leftind = find(pair(blkind,1)==j);
                    switchNum = floor((leftNum-rightNum)/2);
                    ind_rand = randsample(1:length(leftind),switchNum);
                    for k = ind_rand
                        pair((i-1)*blocklength+leftind(k),[1 2]) = pair((i-1)*blocklength+leftind(k),[2 1]);
                    end
                elseif leftNum < rightNum
                    rightind = find(pair(blkind,2)==j);
                    switchNum = floor((rightNum-leftNum)/2);
                    ind_rand = randsample(1:length(rightind),switchNum);
                    for k = ind_rand
                        pair((i-1)*blocklength+rightind(k),[1 2]) = pair((i-1)*blocklength+rightind(k),[2 1]);
                    end
                end
            end
        end
        
        countL = tabulate(pair(blkind,1));
        countR = tabulate(pair(blkind,2));
        countLR(:,:,i) = [countL(:,2) countR(:,2)];
        time = toc
        if time/60 > 5
            if isempty(find(countL(:,2)==0))==1 & isempty(find(countR(:,2)==0))==1 % more lenient criterion
                start = 0;
            else
                pair(blkind,:) = pair_old;
%                 display(i)
%                 display('re2')
                stnd = 're2';
            end
        else
            if isempty(find(abs(countL(:,2)-countR(:,2))>1)) == 1 
                start = 0;
            else
                pair(blkind,:) = pair_old;
%                 display(i)
%                 display('re1')
                stnd = 're1';
            end
        end
    end
    stndmat{i,1} = stnd;
    
    disp(['Done - Block ' num2str(i) ' out of ' num2str(blockNum)])
end

% Check if back to back presentations
row=[];
pair_ori = pair;
for i = 1:blockNum
    %blkind = (i-1)*blocklength+1:i*blocklength;
    c = 1;
    while c~= 0
        c = 0;
        for k = 2:blocklength %length(blkind)
            ind = k+(i-1)*blocklength;
            curr = pair(ind,:);
            pre = pair(ind-1,:);
            if k ~= blocklength                
                if nnz(curr==pre) ~= 0 | nnz(curr==pre([2 1])) ~= 0 % if back to back presentation
                    pair([ind ind+1], :) = pair([ind+1 ind], :); % current one and next one switch
                    c = c + 1;
                    row=[row; ind];
                end
            else
                if nnz(curr==pre) ~= 0 | nnz(curr==pre([2 1])) ~= 0
                    pair([ind ind-blocklength+1], :) = pair([ind-blocklength+1 ind], :);
                    c = c + 1;
                    row=[row; ind];
                end
            end
        end
    end
end

if ~exist(output_dir,'dir') 
    mkdir(output_dir);
end
cd(output_dir) 
save(['tp_' subID '.mat'],'pair')
disp(['Subject ' subID ' done.'])
cd(root_dir) 

end
