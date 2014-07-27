function [] = h__allfiles_collect_features()

yin_path        ='C:\Users\yonatan\Documents\ucsd\tools\yin_changed\';
data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
out_supdir      = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa_processing';

listfile        = [data_supdir,filesep,'uiowa.csv'];
segments_file   = [data_supdir,filesep,'segmentations.mat'];

[metadata,colinds]  = read_list_file(listfile);
load(segments_file);

N                   = size(metadata,1);
capacity            = 4000;
feature             = 'pitch';

switch feature
    case 'pitch'
        feat_dim            = 10;
        outfile             = [out_supdir,filesep,'pitch.mat'];
        indir               = [out_supdir,filesep,'pitch'];
    case 'temporal'
        feat_dim            = 32;
        outfile             = [out_supdir,filesep,'temporal_gainfb.mat'];
        indir               = [out_supdir,filesep,'temporal'];
end

collection          = zeros(feat_dim,capacity);
file_inds           = zeros(1,capacity);
seg_inds            = zeros(1,capacity);

count               = 0;
for ii = 1:N
    filename        = metadata{ii,colinds.filename};
    [path,filecore,ext] = fileparts(filename);
    
    segmentation    = segmentations{ii};
    if ~isstruct(segmentation)
        continue;
    end

    fprintf('%d) %s ------------------\n',ii,filecore);

    class           = metadata{ii,colinds.class};
    type            = metadata{ii,colinds.type};
    out_core        = [class,'.',type,'.',filecore];
    
    % Go over the segments in this file:
    segments        = segmentation.segments;
    for si = 1:size(segments,1)
        switch feature
            case 'pitch'
                infile      = [indir,filesep,out_core,'.seg' num2str(si),'.pitch.mat'];
            case 'temporal'
                infile      = [indir,filesep,out_core,'.seg' num2str(si),'.temporal_gainfb.mat'];
        end
        
        if exist(infile,'file')
            load(infile);
            switch feature
                case 'pitch'
                    to_add  = pitch_features;
                case 'temporal'
                    to_add  = temporal_features;
            end
            fprintf('seg%d,',si);
        else
            fprintf('-- skipping missing file: %s\n',infile);
            continue;
        end
        
        % Store these new filters:
        add         = size(to_add,2);
        i1          = count+1;
        i2          = count+add;
        count       = i2;
        % Do we need more memory allocated?
        if i2 > size(collection,2)
            aloc    = 1000;
            fprintf('== allocate %d more frames\n',aloc);
            collection  = [collection,zeros(size(collection,1),aloc)];
            file_inds   = [file_inds,zeros(1,aloc)];
            seg_inds    = [seg_inds,zeros(1,aloc)];
        end
        collection(:,i1:i2) = to_add;
        file_inds(i1:i2)    = ii;
        seg_inds(i1:i2)     = si;
    end
    
    disp(' ');
    fprintf('count (capacity) = %d (%d)\n',count,size(collection,2));
end

collection  = collection(:,1:count);
file_inds   = file_inds(1:count);
seg_inds    = seg_inds(1:count);
save(outfile,'collection','file_inds','seg_inds');

end

function [representative,gains] = collect_representative_from_segment(seg_file,with_lpc)

load(seg_file);

gains               = -1;

if with_lpc
    % Filter out the low-gain frames:
    [max_gain,ind]  = max(g);
    representative  = A(:,ind);
    gains           = g(ind);
%     min_gain        = min(g);
%     t_gain          = (max_gain + min_gain) / 2;
%     %    q_gain          = quantile(res.g,0.7);
%     frame_inds      = find(g > t_gain);
% 
%     % Choose random example frames:
%     n_take          = min([length(frame_inds),n_ex]);
%     rand_perm       = randperm(length(frame_inds));
%     chosen_fr_inds  = frame_inds(rand_perm(1:n_take));
%     segment_As      = res.A(:,chosen_fr_inds);
else
    powers          = sum(abs(S).^2);
    [max_power,ind] = max(powers);
    representative  = S(:,ind);
end

end
