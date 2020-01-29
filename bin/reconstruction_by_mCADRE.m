function out = reconstruction_by_mCADRE(consmodel,expr_matrix_file,cutoff)

%load('lungExpr.mat');
%load('consRecon2_2.mat');
%cutoff = 75; % resm cutoff
%%
%consmodel = model;
%expr_matrix_file = './data/TCGA_lung.txt';
%cutoff = 75;
%%
expr_pre = file2cell(expr_matrix_file,'\t');
genes = expr_pre(2:end,1);
expr = cell2float(expr_pre(2:end,2:end));
sampleNames = expr_pre(1,2:end);
modela = consmodel;%consistant Recon2
modelgenes = modela.genes;
%genesEntrez = regexp(genes,'([^\.]+)\.\d','tokens');
geneEntrez = modelgenes;
for i  = 1:length(modelgenes)
    tmp = regexp(modelgenes{i},'([^\.]+)\.{0,1}\d{0,1}','tokens');
    geneEntrez{i} = tmp{1}{1};
end
%genesEntrez = [genesEntrez{:}];
%genesEntrez = [genesEntrez{:}];
essentail_rxns = {'biomass_reaction','DM_atp_c_'};
ids = findRxnIDs(modela,essentail_rxns);
out = zeros(length(modela.rxns),length(sampleNames));
cc = 1:length(modela.rxns);
%%
for i  = 1:size(expr,2)
    is_expr = expr(:,i)>cutoff;
    coreGenes = genes(is_expr);
    ExprState = ismember(geneEntrez,coreGenes);
    [rxnState,core_rxns] = gene_stat2rxns_stat_by_rules(modela,ExprState);
    core_rxns = unique([core_rxns;ids']);
    ubiquityScore = zeros(length(modela.rxns),1);
    ubiquityScore(core_rxns) = 10;
    confscore = zeros(length(modela.rxns),1);
    tmodel = mCADRE(modela, ubiquityScore,confscore);
    out_rxns = find(ismember(modela.rxns,tmodel.rxns));
    out(:,i) = ismember(cc,out_rxns);
end
header = ['rxns',sampleNames];
writetxt([header;[modela.rxns,num2cellstr(out)]],'./out/reconstructed_model_mCADRE.txt','\t');