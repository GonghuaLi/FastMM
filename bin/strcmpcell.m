function L = strcmpcell(c,str)
%Return logic for cell _strcmp
L = zeros(1,length(c));
for i  =  1:length(c)
    if strcmp(c{i},str)
        L(i)=1;
    end
end
L=logical(L);
        