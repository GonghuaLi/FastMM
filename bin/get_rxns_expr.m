function rxnState = get_rxns_expr(modelin,enzyme_expr)
% reconstruction Cb model through a given expressed genelist.
% by Gonghua Li
% 2014-04-17
%modelin = liver;
%genstat = gene_stat;
%genelist = modelin.genes(genstat>0);
genes = modelin.genes;
%genesEntrez = regexp(genes,'([^\.]+)\.\d','tokens');
%genesEntrez = [genesEntrez{:}];
%genesEntrez = [genesEntrez{:}];
grRules = modelin.grRules;
%%ExprState = ismember(genesEntrez,genelist);
ExprState = enzyme_expr;%%ismember(genes,genelist);
rxnState  = zeros(length(modelin.rxns),1);
%%
for j = 1:length(grRules)
    if isempty(grRules{j})
        rxnState(j) = -1000;
    else
        line = grRules{j};
        if line(1) ~= '(';
            idxgene = strcmp(line,genes);
            rxnState(j) = ExprState(idxgene);
        else
            positionBrackets = regexp(line,'[\(|\)]');
            score0 = [];
            for k = 1:2:length(positionBrackets)
                charInBrackets = line(positionBrackets(k):positionBrackets(k+1));
                cc = regexp(charInBrackets,'[\w\-]+\.\d','match');
                idxgene = ismember(genes,cc);
                vector_expr = ExprState(idxgene);
                vector_expr(vector_expr<0)=[];
                if regexp(charInBrackets,'and')
                    score = min(vector_expr);
                else
                    if isempty(vector_expr)
                        score = [];
                    else
                        score = sum(vector_expr);
                    end
                end
                if k==1
                    score0 = score;
                elseif positionBrackets(k)- positionBrackets(k-1) == 6
                    if ~isempty(score) & ~isempty(score0) 
                        score0 = min(score0,score);
                    elseif ~isempty(score) & isempty(score0)
                        score0 = score;
                    end                      
                else
                    if ~isempty(score) & ~isempty(score0) 
                        score0 = score0+score;
                    elseif ~isempty(score) & isempty(score0)
                        score0 = score;
                    end  
                end
            end
            if isempty(score0)
                rxnState(j) = -1000;
            else
                rxnState(j) = score0;
            end
        end
    end
end
%core_rxns = find(rxnState>0.01);
%%rxnRemoves = modelin.rxns(rxnState<0.01);
%%modelout = changeRxnBounds(modelin,rxnRemoves,0,'b');
%%modelout.geneStates = ExprState;