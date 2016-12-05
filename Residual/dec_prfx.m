function [level_prfx,i]= dec_prfx(bits,i)

level_prfx = 0;

while i<=length(bits)
    switch bits(i)
        case '0'
            level_prfx = level_prfx + 1;
            i = i + 1;
        case '1'
            i = i + 1;
            return;
    end        
end
