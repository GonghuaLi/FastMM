function flux = FastMM_doubleMetKO(model)
t = clock;
c = ['doubleMetKO',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
cout = [c,'.txt'];
a = cobra2FastKO(model,c);
system(['doubleMetKO -m ',c,' -t max -o ',cout]);
flux = file2cell(cout,'\t');
flux = cell2float(flux);

delete(cout);
rmdir(c,'s');