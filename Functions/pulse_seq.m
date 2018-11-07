function [P] = pulse_seq(x)
% Return structure of all the pulses in a sequence of arbitrary height

% Arbitrary threshold
x = x(:);
thresh = max(x)/10;

pulses = x >= thresh;
trans = diff(pulses);
up = find(trans == 1)+1;
down = find(trans == 1)+1;

if isempty(up) || isempty(down)
    P = [];
else
    range = [up,down];
    dur = (down-up);
    isi = [down(1:end-1),up(2:end)];

    P = struct_from_list('range', range, 'dur', dur, 'isi', isi);
end


