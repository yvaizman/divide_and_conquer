function [] = h__fix_files_wav_format()

data_supdir         = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
listfile            = [data_supdir,filesep,'uiowa.csv'];

[metadata,colinds]  = read_list_file(listfile);

N                   = size(metadata,1);
for ii = 1:N
    wav_file        = [data_supdir,filesep,metadata{ii,colinds.filename}];
    try
        w           = wavread(wav_file);
    catch ME
        sprintf('%d) Fixing file %s',ii,wav_file);
        fix_wav_format(wav_file)
    end
end

end

function [metadata,colinds] = read_list_file(listfile)

[a,b,metadata]  = xlsread(listfile);
headlines       = metadata(1,:);
metadata        = metadata(2:end,:);

colinds         = struct();

for colind = 1:size(headlines,2)
    colinds.(headlines{colind})     = colind;
end

end

function [] = fix_wav_format(wav_file)
fid = fopen(wav_file,'r+');
fseek(fid,20,0);
fwrite(fid,[3 0]);
fclose(fid);
end
