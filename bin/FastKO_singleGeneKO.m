function flux = FastKO_singleGeneKO(model)
t = clock;
c = ['singleGeneKO',num2str(ceil(rand(1,1)*1000)),num2str(ceil(t(6)*100))];
cout = [c,'.txt'];
a = cobra2FastKO(model,c);
system(['singleGeneKO -m ',c,' -t max -o ',cout]);
flux = file2cell(cout,'\t');
flux = cell2float(flux);

delete(cout);
rmdir(c,'s');