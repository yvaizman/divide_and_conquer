function [] = h__allfiles_calculate_spectra_and_filters()

yin_path        = 'C:\Users\yonatan\Documents\ucsd\tools\yin_changed\';
data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
out_supdir      = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa_processing';

listfile        = [data_supdir,filesep,'uiowa.csv'];
segments_file   = [data_supdir,filesep,'segmentations.mat'];

[metadata,colinds]  = read_list_file(listfile);
load(segments_file);

M                   = 14;
params.sr           = 22050;
params.M            = M;
params.hoplen       = 1024;
params.winlen       = 2*params.hoplen;
params.do_preemph   = false;
params.do_yin       = false;
params.yin_path     = yin_path;
params.use_window   = true;

N                   = length(segmentations);

no_lpc_dir          = [out_supdir,filesep,'nolpc'];
lpc_l2_dir          = [out_supdir,filesep,'lpcl2'];
lpc_l1_dir          = [out_supdir,filesep,'lpcl1'];
if ~params.do_preemph
    lpc_l2_dir      = [lpc_l2_dir, '_nopreemph'];
    lpc_l1_dir      = [lpc_l1_dir, '_nopreemph'];
end

if ~exist(no_lpc_dir,'dir')
    mkdir(no_lpc_dir);
end
if ~exist(lpc_l2_dir,'dir')
    mkdir(lpc_l2_dir);
end
if ~exist(lpc_l1_dir,'dir')
    mkdir(lpc_l1_dir);
end

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
    
    wav_file        = [data_supdir,filesep,filename];
    [w,sr_orig]     = wavread(wav_file);
    w               = mean(w,2);
    w               = resample(w,params.sr,sr_orig);

    % Go over the segments in this file:
    segments        = segmentation.segments;
    for si = 1:size(segments,1)
        out_nolpc   = [no_lpc_dir,filesep,out_core,'.seg' num2str(si),'.nolpc.spec.mat'];
        if params.do_preemph
            out_lpcl2   = [lpc_l2_dir,filesep,out_core,'.seg' num2str(si),'.lpcl2M' num2str(M) '.filters.mat'];
            out_lpcl1   = [lpc_l1_dir,filesep,out_core,'.seg' num2str(si),'.lpcl1M' num2str(M) '.filters.mat'];
        else
            out_lpcl2   = [lpc_l2_dir,filesep,out_core,'.seg' num2str(si),'.lpcl2M' num2str(M) '_nopreemph.filters.mat'];
            out_lpcl1   = [lpc_l1_dir,filesep,out_core,'.seg' num2str(si),'.lpcl1M' num2str(M) '_nopreemph.filters.mat'];
        end
        
        seg         = w(segments(si,1):segments(si,2));
        
        if exist(out_nolpc,'file')
            fprintf('-- skipping existing file: %s\n',out_nolpc);
        else
            S       = spectrogram(seg,params.winlen,params.hoplen,params.winlen,params.sr);
            save(out_nolpc,'S');
            fprintf('++ saved no-lpc file: %s\n',out_nolpc);
        end
        
        if exist(out_lpcl2,'file')
            fprintf('-- skipping existing file: %s\n',out_lpcl2);
        else
            params.criterion    = 'L2';
            res     = lpc_analysis_by_frames(seg,params);
            A       = res.A;
            g       = res.g;
            save(out_lpcl2,'A','g');
            fprintf('++ saved lpc-l2 file: %s\n',out_lpcl2);
        end
        
        if exist(out_lpcl1,'file')
            fprintf('-- skipping existing file: %s\n',out_lpcl1);
        else
            params.criterion    = 'L1';
            res     = lpc_analysis_by_frames(seg,params);
            A       = res.A;
            g       = res.g;
            save(out_lpcl1,'A','g');
            fprintf('++ saved lpc-l1 file: %s\n',out_lpcl1);
        end
        
    end
    
end

end

