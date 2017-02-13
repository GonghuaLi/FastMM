function flux = FastMM_doubleGeneKO(model,varargin)
t = clock;
c = ['doubleGeneKO',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
cout = [c,'.txt'];
a = cobra2FastKO(model,c);
isobj = 0;
iscons = 0;

programm_name = 'doubleGeneKO';

if isempty(varargin)
    system([programm_name,' -m ',c,' -t max -o ',cout]);
else
    addition = ' ';
    obj = regexp(varargin{1},'-f\s+([^\s]+)','tokens');
    if ~isempty(obj)
        obj = obj{1};
        obj = obj{1};
        ot = file2cell(obj,'\t');
        ids = findRxnIDs(model,ot(:,1));
        objf = ['obj',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
        writetxt(num2cellstr(ids),objf,'\t');
        addition = [addition,' -f ',objf];
        isobj = 1;
    end
    constrain = regexp(varargin{1},'-c\s+([^\s]+)','tokens');
    if ~isempty(constrain)
        constrain = constrain{1};
        constrain = constrain{1};
        ot = file2cell(constrain,'\t');
        ids = findRxnIDs(model,ot(:,1));
        cons = ['cons',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
        writetxt([num2cellstr(ids),ot(:,2)],cons,'\t');
        addition = [addition,' -c ',cons];
        iscons = 1;
    end   
    sens = regexp(varargin{1},'(-t\s+[^\s]+)','tokens');
    if ~isempty(sens)
        sens = sens{1};
        sens = sens{1};
        addition = [addition,' ',sens];
        system([programm_name,' -m ',c,' -o ',cout,addition]);
    else
        system([programm_name,' -m ',c,' -t max -o ',cout,addition]);
    end
end
if isobj==1
    delete(objf);
end
if iscons==1
    delete(cons);
end
flux = file2cell(cout,'\t');
flux = cell2float(flux);

delete(cout);
rmdir(c,'s');