function [rxnState,core_rxns] = gene_stat2rxns_stat_by_rules(modelin,genstat)
%%
genelist = modelin.genes(genstat>0);
genes = modelin.genes;
%genesEntrez = regexp(genes,'([^\.]+)\.\d','tokens');
%genesEntrez = [genesEntrez{:}];
%genesEntrez = [genesEntrez{:}];
rules = modelin.rules;
%%ExprState = ismember(genesEntrez,genelist);
ExprState = ismember(genes,genelist);
rxnState  = zeros(length(modelin.rxns),1);
x = ExprState;
for i  = 1:length(rxnState)
    if isempty(rules{i})
        rxnState(i) = nan;
    else
        rxnState(i) = eval(rules{i});
    end
end
core_rxns = find(rxnState>0.01);