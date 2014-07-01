function [] = h__mark_separate_notes()

yin_path        = 'C:\Users\yonatan\Documents\ucsd\tools\yin_changed_64\';
addpath(genpath(yin_path));

sr                  = 22050;

yin_params.sr       = sr;
yin_params.wsize    = 2048;
yin_params.hop      = 1024;
smoother            = hamming(10);
smoother            = smoother / mean(smoother.^2);

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
segments_file   = [data_supdir,filesep,'segmentations.mat'];
listfile        = [data_supdir,filesep,'uiowa.csv'];

[metadata,colinds]  = read_list_file(listfile);

N                   = size(metadata,1);
if exist(segments_file,'file')
    load(segments_file);
else
    segmentations   = cell(N,1);
    for ii = 1:N
        segmentations{ii}   = NaN;
    end
    is_bad_file     = false(N,1);
%    wrong_title     = false(N,1);
    fixed_notes     = cell(N,1);
    filenames       = metadata(:,colinds.filename);
end


for ii = 1:N    
    filename        = filenames{ii};
    if any(strcmp(metadata{ii,colinds.class},{'percussion','misc'}))
        % Skip for now
        fprintf('%d) -- skipping perc/misc file: %s\n',ii,filename);
        continue;
    end
    if any(strcmp(metadata{ii,colinds.type},'cello'))%{'viola','viola2012'}))%,'cello'}))
        % Skip for now
        fprintf('%d) -- skipping viola file: %s\n',ii,filename);
        continue;
    end
    if isstruct(segmentations{ii})
        % Then this file was already segmented
        fprintf('%d) -- skipping segmented file: %s\n',ii,filename);
        continue;
    end
    if is_bad_file(ii)
        % Then don't try this file again:
        fprintf('%d) -- skipping bad file: %s\n',ii,filename);
        continue;
    end

    wav_file        = [data_supdir,filesep,filename];
    
    [w,sr_orig]     = wavread(wav_file);
    w               = mean(w,2);
    w               = resample(w,sr,sr_orig);
    L               = length(w);

    % Expected note frequencies:
    if ~isempty(fixed_notes{ii})
        notes       = fixed_notes{ii};
    else
        notes       = metadata{ii,colinds.notes};
    end
    fprintf('%s|%s\n',filename,notes);
    instrument      = metadata{ii,colinds.type};
    [exp_fs,exp_relocts]    = get_note_frequencies(notes,instrument);
    
    % This file might not be tonal, or might only have a single note.
    % In both cases, there is only one segment:
    if length(exp_relocts) <= 1
        segmentations{ii}           = struct();
        segmentations{ii}.segments  = [1,L];
        segmentations{ii}.fs        = exp_fs;
        segmentations{ii}.relocts   = exp_relocts;
        fprintf('%d) ++ single segment file: %s\n',ii,filename);
        save(segments_file,'segmentations','filenames','is_bad_file','fixed_notes');
        continue;
    end
    
    clear segments;
    try
        [segments]                  = segment_to_notes(w,yin_params,smoother,exp_relocts);
        if length(segments) <= 0
            error('bad segmentation');
        end
        segmentations{ii}           = struct();
        segmentations{ii}.segments  = segments;
        segmentations{ii}.fs        = exp_fs;
        segmentations{ii}.relocts   = exp_relocts;
        is_bad_file(ii)             = false;
        fprintf('%d) ++ multi-segment file: %s\n',ii,filename);
    catch ME
        is_bad_file(ii)             = true;
        fprintf('%d) !!! trouble file: %s\n',ii,filename);
    end
    
    % Save the segmentations:
    save(segments_file,'segmentations','filenames','is_bad_file','fixed_notes');
    
end


end

function [segments] = segment_to_notes(w,yin_params,smoother,exp_relocts)

L               = length(w);
% Estimate f0 and powers:
yinr            = yin(w,yin_params);
% Smooth powers:
powers          = filter(smoother,1,yinr.pwr);
est_relocts     = yinr.f0; %YIN provides f0 in units of relative octave
[peaks,locs]    = findpeaks(powers);
peak_relocts    = est_relocts(locs);

% Find representative locations for each of the expected notes:
n_notes         = length(exp_relocts);
center_locs     = zeros(1,n_notes);
loci            = 0;
for ni = 1:n_notes
    reloct      = exp_relocts(ni);
    diffs       = peak_relocts - reloct;
    diffs(1:loci)       = Inf; % Look forward than the latest location
    [mindiff,loci]      = min(abs(diffs));
    if mindiff>(1/14)
        n_notes = 0;
        break;
%        error('didnt find close match for expected frequency');
    end
%    loci        = find(abs(diffs)<(1/25),1,'first');
    center_locs(ni)     = locs(loci);
    % Sweep through the next frames to conquere the whole note area:
    while(loci<length(locs))
        this_loc        = locs(loci);
        next_loci       = loci+1;
        next_loc        = locs(next_loci);
        relocts_between = est_relocts(this_loc:next_loc);
        diffs_between   = relocts_between - reloct;
        % If between the two peaks there is any point with estimated
        % relative octave which is extremely far from the expected, we'll
        % consider the current location the last peak of this note.
        % Otherwise, we'll include the next peak location in the same note:
        if max(abs(diffs_between)) < 1/30
            loci        = next_loci;
            % Mark the highest-power point with the same note as the
            % "center" location:
            if powers(next_loc)>powers(center_locs(ni))
                center_locs(ni) = next_loc;
            end
        else
            break;
        end
    end
end

center_samples  = yin_params.hop*(center_locs-1) + 1;
squares         = w.^2;
sample_powers   = filter(ones(50,1),1,squares);
segments        = zeros(n_notes,2);
last_end        = 0;
% Find power valeys between notes:
for ni = 1:n_notes
    from            = last_end + 1;
    segments(ni,1)  = from;
    if ni == n_notes
        segments(ni,2)  = L;
        break;
    end
    
    this_center     = center_samples(ni);
    next_center     = center_samples(ni+1);
    search_powers   = sample_powers;
    search_powers([1:this_center,next_center:end])  = Inf;
    [valey,sample]  = min(search_powers);
    last_end        = sample; %(loc-1)*yin_params.hop;
    segments(ni,2)  = last_end;
end

%helper_plot2(w,segments);
end

function [uniform_segments] = get_uniform_segments(n_notes,L)
tt = round((0:(n_notes-1))*L/n_notes)'+1;
uniform_segments = [tt,[tt(2:end)-1;L]];
end

function [] = helper_plot2(w,segments)

figure;
plot(w,'b');
minval = min(w);
maxval = max(w);
for si = 1:size(segments,1)
    hold on;
    plot([segments(si,1),segments(si,1)],[minval,maxval],'.-g');
    hold on;
    plot([segments(si,2),segments(si,2)],[minval,maxval],'o-k');
end

end

function [] = helper_plot(est_relocts,center_locs,exp_relocts,powers)

figure;
plot(est_relocts,'b');
hold on;
plot(center_locs,exp_relocts,'og');
hold on;plot(4+powers*400,'m')

end

function [fs,relocts] = get_note_frequencies(notes,instrument)

if isnan(notes)
    fs              = [];
    relocts         = [];
    return;
end

parts               = regexp(notes,'_','split');
notes               = parts{1};
extra_notes         = parts(2:end);

if length(notes) < 4
    % Then there is just one note.
    [fs,relocts]    = note2freq(notes);
    return;
end

% Look for the number of the first note:
is_digit            = isdigit(notes);
ind                 = find(is_digit,1);

note1               = notes(1:ind);
note2               = notes((ind+1):end);

[f1,reloct1]        = note2freq(note1);
[f2,reloct2]        = note2freq(note2);

relocts             = sort(union(reloct1:(1/12):reloct2,reloct2),'ascend');

% Handle extra notes (probably quarter-tones):
if ~isempty(extra_notes)
    for ei = 1:length(extra_notes)
        after_ind       = str2num(extra_notes{ei});
        extra_reloct    = mean([relocts(after_ind),relocts(after_ind+1)]);
        relocts         = [relocts,extra_reloct];
    end
    relocts             = sort(relocts,'ascend');
end

if strcmp('viola2012',instrument)
    relocts             = relocts - (1+(1/12));
end

fs                  = 440 * 2.^relocts;

end

function y = isdigit(x)
%ISDIGIT	True for digit characters.
%
% y = isdigit(x)
%    Given a string x, return a value y the same shape as x with 1's where x
%    has digit characters and 0's elsewhere.  Digits are simply '0' to '9'.
%    x need not be a string, merely have values in the right range.
%
% See also ISLETTER, ISSPACE, ISSTR, IS2POWER, ISEMPTY, 

y = (x >= '0') & (x <= '9');

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
