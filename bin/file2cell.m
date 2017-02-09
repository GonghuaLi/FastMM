function lines = file2cell(filename,varargin)
%%%%%%
%usage c = file2cell(filename,pat)
%      c = file2cell(filename)
%Read file data to matlab cell arrays.
%Writen by Gong-Hua Li
%2009-09-29
%
txt = fileread(filename);
%del '\r' in windows system
%ridx = findstr(txt,sprintf('\r'));
%txt(ridx)=[];
txt = strrep(txt,sprintf('\r'),'');
lines = regexp(txt,'\n+','split');
clear txt;
if strcmp(lines{end},'')
    lines(end)=[];
end
if isempty(lines)
    return;
end

if isempty(varargin) %read file to 1xn cell array
    return;
elseif length(varargin)==1 %read file to nxm cell array
    len = length(lines);
    p=regexp(lines,varargin{1},'split');
    maxpat = length(p{1});
    minpat = maxpat;
    for i  = 1:length(p)
        maxpat = max(length(p{i}),maxpat);
        minpat = min(length(p{i}),minpat);
    end
    if minpat == maxpat
        lines = reshape([p{:}],maxpat,len);
        lines = lines';
    else
        warning(['Eachline elements is not the same, automatedly add: ''',varargin{1},'''  in the end.']); %#ok<WNTAG>
        lines = repmat({''},length(lines),maxpat);
        for i = 1:length(p)
             lines(i,1:length(p{i}))=p{i};
        end
    end
else
    error('Erro:The input must be one or two argvs!');
end