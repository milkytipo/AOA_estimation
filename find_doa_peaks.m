function doa = find_doa_peaks(spec,theta,N)
    
    [pks,loc] = findpeaks(spec);
    [V,I] = sort(pks,'descend');
    if ~isempty(pks)
        doa = sort(theta(loc(I(1:N))),'ascend');
    else
        [mv,mi] = max(spec);
        doa = theta(mi);
    end
end