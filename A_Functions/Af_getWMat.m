function W = Af_getWMat(filename)
% This function will get the W matrix from a .hh wind file
%
% Input:        filename - name of file

%% Find First line of input
fid = fopen([filename,'.wnd']);
line = '!'; iLine = 0;
while strcmp(line(1),'!')
    iLine = iLine + 1;
    line = fgetl(fid);
end
fclose(fid);

filename    = [filename,'.wnd'];             % append .hh (for now)
W           = dlmread(filename,'\t',iLine-1,0);
