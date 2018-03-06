if ~isempty(strfind(param.hostname, 'Justins-Mac')) | ...
   ~isempty(strfind(param.hostname, 'willis'))
	p = gcp('nocreate');
else
	p = gcp;
end
