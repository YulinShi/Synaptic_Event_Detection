function TE=center_mass(start_pos, end_pos, fitting)
up=0;
down=0;
for i= start_pos:end_pos
    up=up+i*fitting(i);
    down=down+fitting(i);
end
TE=round(up/down);