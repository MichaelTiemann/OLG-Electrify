function F=Electrify_GE_wages(L_h,L_f,w)
% A continuous function that returns zero when L_h is balanced with L_f
% and otherwise attempts to steer wages (w) toward balancing L_h and L_f

F=0;

if (L_h-L_f)^2 <= 10^(-6)
    return
elseif L_f>L_h
    % If w>1, this attenuates overdemand by firms
    % If w<1, this amplifies overdemand by firms
    F=(L_f/L_h)/w;
else
    % If w>1, this reverses and amplifies oversupply by households
    % If w<1, this reverses but attenuates oversupply by households
    F=-(L_h/L_f)*w;
end

end