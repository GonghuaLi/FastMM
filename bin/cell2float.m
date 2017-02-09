function s = cell2float(c)
%Convert cell array to numeric array
c(strcmp('',c))={'0'};
space = ' ';
p=char(c);
p = [p,space(ones(size(p,1),1))];
p = p';
s = str2float(p(:));
s = reshape(s,size(c));