function [f,oct_relA4] = note2freq(note,instrument)

octave      = str2num(note(end));
oct_rel4    = octave - 4;
letter      = note(1:end-1);
relA        = relation_to_A_within_octave(letter);

oct_relA4   = oct_rel4 + relA/12;

f           = 440 * 2^oct_relA4;

end

function [relA] = relation_to_A_within_octave(letter)

map.C   = -9;
map.Db  = -8;
map.D   = -7;
map.Eb  = -6;
map.E   = -5;
map.F   = -4;
map.Gb  = -3;
map.G   = -2;
map.Ab  = -1;
map.A   = 0;
map.Bb  = 1;
map.B   = 2;

relA    = map.(letter);

end