function [] = h__testing()

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\wav';

wav_file        = [data_supdir,filesep,'piano\piano\Piano.ff.Eb4.wav'];

sr              = 22050;
M               = 20;
hoplen          = 1024;
winlen          = 2*hoplen;
do_emph         = true;

[w,sr_orig]     = wavread(wav_file);
w               = mean(w,2);
w               = resample(w,sr,sr_orig);

[A,g,e,preemph] = lpc_analysis_by_frames(w,sr,M,do_emph,hoplen,winlen);

end