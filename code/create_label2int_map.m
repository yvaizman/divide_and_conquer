function [label2int,int2label] = create_label2int_map(labels)

int2label   = cell(1,0);
label2int   = struct();
for li = 1:length(labels)
    label   = labels{li};
    if ~any(strcmp(label,int2label))
        int2label{1,end+1}  = label;
        label2int.(label)   = length(int2label);
    end
end

end