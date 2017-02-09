function a = cobra2FastKO(inmodel,outmodelname)
%% contain:
% .mps  file: linner programm model for cflux
% .rxns file: reaction annotation file for cflux
% .mets file:  mets annotation file for cflux
% .gens  file:  gene information file for cflux
% .ruls file: gene rules for reactions
% .objf  file: objective functions file for cflux

if ~isdir(outmodelname)
    mkdir(outmodelname)
end
if regexp(outmodelname,'\/')
    sp = '/';
    filename = regexp(outmodelname,'\/','split');
    if isempty(filename{end})
        filename = filename{end-1};
    else
        filename = filename{end};
    end
elseif regexp(outmodelname,'\\')
    sp = '\';
    filename = regexp(outmodelname,'\\','split');
    if isempty(filename{end})
        filename = filename{end-1};
    else
        filename = filename{end};
    end
else
    sp = '/';
    filename = outmodelname;
end
%% .mps
a = cobra2mps(inmodel,[outmodelname,sp,filename,'.mps']);

%% .rxns
fid = fopen([outmodelname,sp,filename,'.rxns'],'wb');
for i = 1:length(inmodel.rxns)
    fprintf(fid,'%5d %d %s %s %s\n',i,inmodel.rev(i),formatNames('R',i),...
        formatDesnames(inmodel.rxns{i}),inmodel.rxnNames{i});
end
fclose(fid);
%% .mets
fid = fopen([outmodelname,sp,filename,'.mets'],'wb');
for i = 1:length(inmodel.mets)
    fprintf(fid,'%5d %d %s %s %s\n',i,1,formatNames('M',i),...
        formatDesnames(inmodel.mets{i}),inmodel.metNames{i});
end
fclose(fid);

%% .gens
fid = fopen([outmodelname,sp,filename,'.gens'],'wb');
for i = 1:length(inmodel.genes)
    fprintf(fid,'%5d %d %s %s\n',i,1,formatNames('G',i),...
        formatDesnames(inmodel.genes{i}));
end
fclose(fid);

%% .ruls
fid = fopen([outmodelname,sp,filename,'.ruls'],'wb');
for i = 1:length(inmodel.rules)
    outstr = regexprep(inmodel.rules{i},'x\((\d+)\)','$1');
    if length(outstr)>1024
        outstr = '';
    end
    fprintf(fid,'%5d %d %s %s\n',i,1,formatNames('R',i),...
        outstr);
end
fclose(fid);
%% .obj
fid = fopen([outmodelname,sp,filename,'.objf'],'wb');
aa  = find(inmodel.c>0);
for i = 1:length(aa)
    fprintf(fid,'%d\n',aa(i));
end
fclose(fid);




function c = formatNames(rtype, b)
cb = num2str(b);
c = rtype;
for i = 1:5-length(cb)
    c = [c,'0'];
end
c = [c,cb];

function c = formatDesnames(str)
if length(str)>=50
    c = str(1:50);
else
    c = str;
    c = [c,repmat(' ',1,50-length(str))];
end



