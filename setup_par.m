if ~isempty(strfind(param.hostname, 'Justins-Mac')) | ...
   ~isempty(strfind(param.hostname, 'willis'))
	p = gcp('nocreate');
%     	if isempty(p)
%         	parpool(2);
%     	elseif ~p.NumWorkers==2
%         	delete(p);
%         	parpool(2);
%         	p = gcp('nocreate');
%     	end
else
	p = gcp;
end
