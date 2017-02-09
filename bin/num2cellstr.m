function c = num2cellstr(d)
[m,n] = size(d);
out1 = num2str(d)';
out1 = [out1;repmat(' ',1,m)];
out2 = out1(:)';
c = regexp(out2,'\s+','split');
id = ones(1,length(c));
if isempty(c{1})
    id(1) = 0;
end
if isempty(c{end})
    id(end) = 0;
end
c = c(logical(id));
c = reshape(c,n,m)';
