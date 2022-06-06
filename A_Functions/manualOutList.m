function OutList = manualOutList(sumFileName)

if 0
    sumFileName = './FAST8_IF/FASTinput2.SFunc.sum';
end

fid = fopen(sumFileName);

lookingFor =  '   Requested Channels in FAST Output File(s)  ';
OutList = [];
iOutput = 1;


tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    
    if strcmp(tline,lookingFor)         % start of channel list
        for iLine = 1:4                 % collect 4 more lines
            tline = fgetl(fid);
        end
        
        while ischar(tline)
            lineData = textscan(tline,'%d\t%s\t%s\t%s');
            OutList = [OutList;lineData{2}];
            tline = fgetl(fid);
        end
    end
    
end

fclose(fid);

