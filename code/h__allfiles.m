function [] = h__allfiles()

yin_path        ='C:\Users\yonatan\Documents\ucsd\tools\yin_changed\';
data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
listfile        = [data_supdir,filesep,'uiowa.csv'];
segments_file   = [data_supdir,filesep,'segmentations.mat'];

[metadata,colinds]  = read_list_file(listfile);
load(segments_file);

is_percussion       = strcmp('percussion',metadata(:,colinds.class));
is_misc             = strcmp('misc',metadata(:,colinds.class));

params.sr           = 22050;
params.M            = 20;
params.hoplen       = 1024;
params.winlen       = 2*params.hoplen;
params.do_preemph   = true;
params.do_yin       = false;
params.yin_path     = yin_path;
params.use_window   = true;
params.criterion    = 'L2';

N                   = length(segmentations);
capacity            = 10000;
As                  = zeros(params.M+1,capacity);
file_inds           = zeros(1,capacity);

count               = 0;
for ii = 1:N
    filename        = metadata{ii,colinds.filename};
    
    segmentation    = segmentations{ii};
    if ~isstruct(segmentation)
        continue;
    end

    fprintf('%d) %s\n',ii,filename);

    wav_file        = [data_supdir,filesep,filename];
    [w,sr_orig]     = wavread(wav_file);
    w               = mean(w,2);
    w               = resample(w,params.sr,sr_orig);

    % Go over the segments in this file:
    segments        = segmentation.segments;
    for si = 1:size(segments,1)
        seg         = w(segments(si,1):segments(si,2));
        segment_As  = collect_analyses_from_segment(seg,params);
        % Store these new filters:
        add         = size(segment_As,2);
        i1          = count+1;
        i2          = count+add;
        count       = i2;
        % Do we need more memory allocated?
        if i2 > size(As,2)
            aloc    = 1000;
            fprintf('== allocate %d more frames\n',aloc);
            As      = [As,zeros(size(As,1),aloc)];
            file_inds   = [file_inds,zeros(1,aloc)];
        end
        As(:,i1:i2) = segment_As;
        file_inds(i1:i2)    = ii;
    end
    
    fprintf('count (capacity) = %d (%d)\n',count,size(As,2));
end

end

function [segment_As] = collect_analyses_from_segment(seg,params)

n_ex            = 1;
res             = lpc_analysis_by_frames(seg,params);

% Filter out the low-gain frames:
max_gain        = max(res.g);
min_gain        = min(res.g);
t_gain          = (max_gain + min_gain) / 2;
%    q_gain          = quantile(res.g,0.7);
frame_inds      = find(res.g > t_gain);

% Choose random example frames:
n_take          = min([length(frame_inds),n_ex]);
rand_perm       = randperm(length(frame_inds));
chosen_fr_inds  = frame_inds(rand_perm(1:n_take));
segment_As      = res.A(:,chosen_fr_inds);

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
