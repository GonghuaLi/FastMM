function L = isemptycell(c)
%Return logic for cell _isempty
L = zeros(1,length(c));
for i  =  1:length(c)
    if isempty(c{i})
        L(i)=1;
    end
end
L=logical(L);