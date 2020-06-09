function [c, time] = run_hmetis(name,balance,parts)
% Matlab interface for calling hmetis and recovering clustering
% from output files
curdir = pwd;
file = strcat('../../data/',name);

cd('../include/hmetis-1.5-osx-i686')
command = strcat('./shmetis',{' '},file,{' '},num2str(parts),{' '},num2str(balance));

tic;
[status,cmdout] = system(command{1});
time = toc;
cmdout
partfile = strcat(file,'.part.',num2str(parts));
fid = fopen(partfile,'r');
c = fscanf(fid,'%d')+1;
fclose(fid);
rmcmd = strcat('rm',{' '},partfile);    % don't need this file around
system(rmcmd{1});
cd(curdir)
end