function stata = writetxt(s,filename,varargin)
%write cell array or string array or numeric array to a file
%writen by lgh
%2009-10-13
if isnumeric(s)
    s = num2str(s);
end
if min(size(s))==1 || ischar(s)
     s=join(s,'\n');
else
    if isempty(varargin)
        pat = '\t';
    else
        pat = varargin{1};
    end
    p=s(:,1);
    for i = 1:size(s,1)
        p{i}=join(s(i,:),pat);
    end
    s=join(p,'\n');
end
fid = fopen(filename,'w');
fprintf(fid,'%s',s);
fclose(fid);
stata=1;





