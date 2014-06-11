function [synth_w] = synthesis_by_frames(A,g,e,preemph,hoplen,winlen)

L           = length(e);
n_frames    = length(g);

synth_w     = zeros(L,1);
window      = hamming(winlen);
for fi = 1:n_frames
    from    = (fi-1)*hoplen + 1;
    to      = from+winlen-1;
    if to > L
        to      = L;
        window  = hamming(to-from+1);
    end
    
    eframe  = e(from:to);
    a       = A(:,fi);
    gain    = g(fi);
    geframe = gain * eframe;
    sframe  = filter(1,a,geframe);
    wframe  = sframe .* window;
    
    synth_w(from:to)    = synth_w(from:to) + wframe;
end

% De-emphasize the pre-emphasis that was done:
synth_w     = filter(1,preemph,synth_w);

end