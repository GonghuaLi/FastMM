%% install
if ismac
    system('chmod 755 ./bin/linux64/*');
    system('cp ./bin/mac/* ./bin/');
    disp('Use mac binary.');
elseif isunix
    system('chmod 755 ./bin/linux64/*');
    system('cp ./bin/linux64/* ./bin/');
    disp('Use linux64 binary'); 
else
    disp('Use win32/64 binary');
    system('copy .\bin\win32\* .\bin\');
end










    
