function a  = cobra2mps(model,outfile)
a =1;
S = model.S;

fid = fopen(outfile,'w');
fprintf(fid,'NAME          CFlux                                                             \n');

fprintf(fid,'ROWS\n');
if isfield(model,'csense')
    for i = 1:length(model.csense)
        fprintf(fid,' %s  %s  \n',model.csense(i), formatNames('M',i));
    end
else
    for i = 1:length(model.mets)
        fprintf(fid,' E  %s  \n',formatNames('M',i));
    end
end

fprintf(fid,'COLUMNS\n');
for i = 1:length(model.rxns)
    indx = full(S(:,i))~=0;
    reaction_mets = find(full(S(:,i))~=0);
    col_value = full(S(indx,i));
    for j = 1:length(reaction_mets)
        if rem(j,2)==1
            fprintf(fid,'    %s    %s%16.5f   ',formatNames('R',i),...
                formatNames('M',reaction_mets(j)),col_value(j));
        else
            fprintf(fid,'%s%16.5f   \n',formatNames('M',reaction_mets(j)),col_value(j));
        end
        if j == length(reaction_mets) & rem(j,2)==1
            fprintf(fid,'\n');
        end
    end
end

tzero = 0;
fprintf(fid,'RHS\n');
for i = 1:length(model.mets)
    if rem(i,2) == 1
        fprintf(fid,'    .00001    %s%16.5f   ',formatNames('M',i),model.b(i));
    else
        fprintf(fid,'%s%16.5f   \n',formatNames('M',i),model.b(i));
    end
    if i == length(model.mets) & rem(i,2)==1
            fprintf(fid,'\n');
    end
end

fprintf(fid,'BOUNDS\n');
for i = 1:length(model.rxns)
    fprintf(fid,' LO BND001    %s%16.6f   \n',formatNames('R',i),model.lb(i));
    fprintf(fid,' UP BND001    %s%16.6f   \n',formatNames('R',i),model.ub(i));
    %fprintf(fid,' LO BND001    %s%16.3f   \n',formatNames('R',i),model.lb(i));
    %fprintf(fid,' UP BND001    %s%16.3f   \n',formatNames('R',i),model.ub(i));
end
fprintf(fid,'ENDATA\n');

fclose(fid);

function c = formatNames(rtype, b)
cb = num2str(b);
c = rtype;
for i = 1:5-length(cb)
    c = [c,'0'];
end
c = [c,cb];





