function [metadata,colinds] = read_list_file(listfile)

[a,b,metadata]  = xlsread(listfile);
headlines       = metadata(1,:);
metadata        = metadata(2:end,:);

colinds         = struct();

for colind = 1:size(headlines,2)
    colinds.(headlines{colind})     = colind;
end

end
