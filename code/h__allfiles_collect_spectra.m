function [] = h__allfiles_collect_spectra()

yin_path        ='C:\Users\yonatan\Documents\ucsd\tools\yin_changed\';
data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
out_supdir      = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa_processing';

listfile        = [data_supdir,filesep,'uiowa.csv'];
segments_file   = [data_supdir,filesep,'segmentations.mat'];

[metadata,colinds]  = read_list_file(listfile);
load(segments_file);

is_percussion       = strcmp('percussion',metadata(:,colinds.class));
is_misc             = strcmp('misc',metadata(:,colinds.class));

M                   = 14;
lpc_version         = 'l2';

params.sr           = 22050;
params.M            = M;
params.hoplen       = 1024;
params.winlen       = 2*params.hoplen;
params.do_preemph   = false;
params.do_yin       = false;
params.yin_path     = yin_path;
params.use_window   = true;

mfcc_params.fs                  = params.sr;
mfcc_params.fft_size            = params.winlen;
mfcc_params.n_mel_bins          = 128;
mfcc_params.num_ceps_coeffs     = 14;
mfcc_params.use_first_coeff     = false;
[mel_mat,dct_mat,c] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);
spec_size           = size(mel_mat,2);

with_lpc            = true;

N                   = length(segmentations);
capacity            = 3500;
if with_lpc
    collection      = zeros(params.M+1,capacity);
    if params.do_preemph
        outfile         = [out_supdir,filesep,'representatives.As.lpc' lpc_version 'M' num2str(M) '.mat'];
    else
        outfile         = [out_supdir,filesep,'representatives.As.lpc' lpc_version 'M' num2str(M) '_nopreemph.mat'];
    end
else
    collection      = zeros(spec_size,capacity);
    outfile         = [out_supdir,filesep,'representatives.Ss.nolpc.mat'];
end
file_inds           = zeros(1,capacity);
seg_inds            = zeros(1,capacity);
gains               = -ones(1,capacity);

no_lpc_dir          = [out_supdir,filesep,'nolpc'];
lpc_dir             = [out_supdir,filesep,'lpc', lpc_version];
if ~params.do_preemph
    lpc_dir         = [lpc_dir, '_nopreemph'];
end

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
        in_nolpc    = [no_lpc_dir,filesep,out_core,'.seg' num2str(si),'.nolpc.spec.mat'];
        if params.do_preemph
            in_lpc      = [lpc_dir,filesep,out_core,'.seg' num2str(si),'.lpc' lpc_version 'M' num2str(M) '.filters.mat'];
        else
            in_lpc      = [lpc_dir,filesep,out_core,'.seg' num2str(si),'.lpc' lpc_version 'M' num2str(M) '_nopreemph.filters.mat'];
        end
        
        if with_lpc
            infile  = in_lpc;
        else
            infile  = in_nolpc;
        end
        
        if exist(infile,'file')
            [to_add,add_gains]  = collect_representative_from_segment(infile,with_lpc);
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
            gains       = [gains,-ones(1,aloc)];
        end
        collection(:,i1:i2) = to_add;
        file_inds(i1:i2)    = ii;
        seg_inds(i1:i2)     = si;
        gains(i1:i2)        = add_gains;
    end
    
    disp(' ');
    fprintf('count (capacity) = %d (%d)\n',count,size(collection,2));
end

collection  = collection(:,1:count);
file_inds   = file_inds(1:count);
seg_inds    = seg_inds(1:count);
gains       = gains(1:count);
save(outfile,'collection','file_inds','seg_inds','gains');

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
