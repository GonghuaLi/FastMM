function flux = FastKO_FVA(model)
t = clock;
c = ['fva',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
cout = [c,'.txt'];
a = cobra2FastKO(model,c);
system(['FVA -m ',c,' -o ',cout]);
flux = file2cell(cout,'\t');
flux = cell2float(flux);

delete(cout);
rmdir(c,'s');