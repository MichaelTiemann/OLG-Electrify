function [sprime,aprime,s,a]=decode_sa_single(saprime,sa,single_0,single_1)
% Helper function to split shares-and-assets into shares and assets.
% I have no idea why we need to cast so many of these values, but we do.

sprime=single_0;
aprime=saprime;
s=single_0;
a=sa;

if saprime>=single_1
    sprime=floor(saprime);
    aprime=rem(saprime,single_1);
end

if sa>=single_1
    s=floor(sa);
    a=rem(sa,single_1);
end


end