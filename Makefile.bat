# glpk win32
c:\mingw\bin\gcc ./src/singleGeneKO.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/singleGeneKO.exe
c:\mingw\bin\gcc ./src/doubleGeneKO.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/doubleGeneKO.exe
c:\mingw\bin\gcc ./src/FVA.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/FVA.exe
c:\mingw\bin\gcc ./src/singleGeneKI.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/singleGeneKI.exe
c:\mingw\bin\gcc ./src/singleMetKO.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/singleMetKO.exe
c:\mingw\bin\gcc ./src/doubleMetKO.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/doubleMetKO.exe
c:\mingw\bin\gcc ./src/FBA.c C:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/FBA.exe


D:\mingw\bin\gcc ./src/singleGeneKO.c D:/MinGW/gnu/lib/libglpk.a -I./include -O3 -lm -o ./bin/win32/singleGeneKOa.exe


# gurobi win64
gcc ./src/FBA_gurobi.c -m64 -g -o ./bin/win32/FBA_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/FVA_gurobi.c -m64 -g -o ./bin/win32/FVA_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/singleGeneKO_gurobi.c -m64 -g -o ./bin/win32/singleGeneKO_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/doubleGeneKO_gurobi.c -m64 -g -o ./bin/win32/doubleGeneKO_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/singleGeneKI_gurobi.c -m64 -g -o ./bin/win32/singleGeneKI_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/singleMetKO_gurobi.c -m64 -g -o ./bin/win32/singleMetKO_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3
gcc ./src/doubleMetKO_gurobi.c -m64 -g -o ./bin/win32/doubleMetKO_gurobi.exe -IC:\gurobi652\win64\include  -I.\include -LC:\gurobi652\win64\lib -lm -lgurobi65 -lpthread -O3

gcc ./src/FVA_gurobi_mat.c -m64 -g -o ./bin/win32/FVA_gurobi_mat -IC:\gurobi652\win64\include -I.\include -I"C:\Program Files\MATLAB\R2013a\extern\include" -LC:\gurobi652\win64\lib -L"C:\Program Files\MATLAB\R2013a\bin\win64" -lmat -leng -lmx -lgurobi65 -lpthread -lm -O3

# cplex win64
gcc ./src/FVA_cplex_mat.c -m64 -g -o ./bin/win32/FVA_cplex_mat -ID:/programms/IBM/ILOG/CPLEX_Studio126/cplex/include/ -I./include -ID:/programms/Matlab/R2017a/extern/include -LD:/programms/IBM/ILOG/CPLEX_Studio126/cplex/lib/x64_windows_vs2010/stat_mda -LD:/programms/Matlab/R2017a/bin/win64 -lmat -leng -lmx -lcplex1260 -lpthread -lm -O3

gcc ./src/singleGeneKO_cplex.c -m64 -g -o ./bin/win32/singleGeneKO_cplex -ID:/programms/IBM/ILOG/CPLEX_Studio126/cplex/include/ -I./include -ID:/programms/Matlab/R2017a/extern/include -LD:/programms/IBM/ILOG/CPLEX_Studio126/cplex/lib/x64_windows_vs2010/stat_mda -LD:/programms/Matlab/R2017a/bin/win64 -lmat -leng -lmx -lcplex1260 -lpthread -lm -O3

gcc ./src/doubleGeneKO_cplex.c -m64 -g -o ./bin/win32/doubleGeneKO_cplex -ID:/programms/IBM/ILOG/CPLEX_Studio126/cplex/include/ -I./include -ID:/programms/Matlab/R2017a/extern/include -LD:/programms/IBM/ILOG/CPLEX_Studio126/cplex/lib/x64_windows_vs2010/stat_mda -LD:/programms/Matlab/R2017a/bin/win64 -lmat -leng -lmx -lcplex1260 -lpthread -lm -O3

gcc ./src/singleGeneKI_cplex.c -m64 -g -o ./bin/win32/singleGeneKI_cplex -ID:/programms/IBM/ILOG/CPLEX_Studio126/cplex/include/ -I./include -ID:/programms/Matlab/R2017a/extern/include -LD:/programms/IBM/ILOG/CPLEX_Studio126/cplex/lib/x64_windows_vs2010/stat_mda -LD:/programms/Matlab/R2017a/bin/win64 -lmat -leng -lmx -lcplex1260 -lpthread -lm -O3
