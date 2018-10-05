function setup_par
doCreatePar = false;
if doCreatePar
   if ~isempty(strfind(param.hostname, 'Justins-Mac')) | ...
         ~isempty(strfind(param.hostname, 'willis'))
      p = gcp('nocreate');
      %           if isempty(p)
      %               parpool(2);
      %       elseif ~p.NumWorkers==2
      %               delete(p);
      %               parpool(2);
      %               p = gcp('nocreate');
      %       end
   else
      pause ( rand * 5 ); % random pause for up to 5 seconds to avoid filename collisions in parallel logging code.
      p = gcp;
   end
end

