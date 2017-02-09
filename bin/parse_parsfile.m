function pars = parse_parsfile(parsfile)
parsInfor = file2cell(parsfile,'\n');
for i  = 1:length(parsInfor)
    thisline = parsInfor{i}
    if strcmp(thisline(1),'#') || length(thisline)<5
        continue;
    end
    m = regexp(thisline,'''([^'']+)''','tokens');
    m = m{1};
    m = m{1};
    if regexp(thisline,'^consModel')
        pars.consModel = m;
    elseif regexp(thisline,'^exchangeFile')
        pars.exchangeFile = m;
    elseif regexp(thisline,'^expressionFile');
        pars.expressionFile = m;
    elseif regexp(thisline,'^expressionType')
        pars.expressionType = m;
    elseif regexp(thisline,'^sampleInfor')
        pars.sampleInfor = m;
    elseif regexp(thisline,'^cutoff')
        pars.cutoff = str2num(m);
    elseif regexp(thisline,'^cpu')
        pars.numCPU = str2num(m);
    elseif regexp(thisline,'^ismcmc')
        pars.ismcmc = m;
    end
end
        
        