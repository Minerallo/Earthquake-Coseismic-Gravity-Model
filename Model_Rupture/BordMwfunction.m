function [minMw,maxMw]=BordMwfunction(button,longitude,lonlimwest,lonlimest,latitude,latlimdown,latlimup,Moexpo,dateRquake,datein,dateout)
    k=find(longitude>=lonlimwest & longitude<=lonlimest & dateRquake>=datein & dateRquake>=dateout);
    %find all the earthquake in longitude limit, latitude limit and Mo limit and Give the index to keep for the concerned zone
    for T=k
        ind1=find(latitude(T)>=latlimdown & latitude(T)<=latlimup & dateRquake(T)>=datein & dateRquake(T)>=dateout);
    end
    if isempty(ind1)==1
        disp('No earthquake found in this zone');
    elseif isempty(ind1)==0
        disp([num2str(numel(ind1)),' earthquakes found in this zone']);
    end
    % Extracted parameters for earthquakes define
    for i=1:numel(ind1)
        Rt=k(ind1)';
        Mo=Moexpo(Rt);
    end
    Mw = (2/3)*log10(Mo*1e7)-10.7; % do the conversion of each Mo of each earthquake in Mw.
    minMw=min(Mw);
    maxMw=max(Mw);
end
