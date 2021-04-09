% clc
% clear all;
% load('mm10_refseq');
% load noncodeproject_10;
% Candidate_transcript = cell(length(Chr),1);
% Bps = zeros(length(Chr),1);
% % 10k cutoff
% cutoff = 5000;
for i = 1:length(Chr) % for each non-coding RNA
    if Strand{i,1} == '+';
        ix = find(strcmp(ChrName,Chr{i,1})&strcmp(Trend,'+')&position2(i)<CdsStart&strncmp(TranscriptName,'NM_',3));% match chr, and match strand, and the ones are downstream
        if ~isempty(ix)
            tmp = CdsStart(ix)-position2(i);
            try 
                [Bps(i),idx] = min(tmp); % report smallest distance and index
            catch
                1
            end
            Candidate_transcript{i,1} = TranscriptName{ix(idx),1};
        else
            Candidate_transcript{i,1} = 0;
        end
        
    else
        ix = find(strcmp(ChrName,Chr{i,1})&strcmp(Trend,'-')&position1(i)>CdsEnd&strncmp(TranscriptName,'NM_',3));
        if ~isempty(ix)
            tmp = position1(i) - CdsEnd(ix);
            [Bps(i),idx] = min(tmp);
            Candidate_transcript{i,1} = TranscriptName{ix(idx),1};
        end
    end
    
end


% CC = zeros(size(testcoexp1,1),length(proteinc));
% Noncode = testcoexp1(1:9,:);
% proteinc = testcoexp1(10:19,:);
% 
% for i= 1:size(Noncode);
%     for j = 1:size(proteinc);
%         tmp = corrcoef(Noncode(i,:),proteinc(j,:));
%         CC(i,j) =tmp(1,2);
%     end
% end
% 
