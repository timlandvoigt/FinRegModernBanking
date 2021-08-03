function [val,key] = replaceVal( var, optional, set, full, mapDict, optionalDefaults )
	optPattern = '\[(.*?)\]';
	tk = regexp(optional,optPattern,'tokens');
	tk = [tk{:}];
	optElements = optionalDefaults;
	overwriteIdx = ~isempty(tk) & ~strcmp(tk,'');
	optElements( overwriteIdx ) = tk( overwriteIdx );
	if strcmp(set,'param')
		key = ['param - ',var,' - ',optElements{1}];
		val = sprintf( '%0.3g', mapDict(key) );
	else
		key = [set,' - ',optElements{3},' - ',optElements{2},' - ',var,' - ',optElements{1}];
		try
			val = sprintf( ['%.',optElements{4},'f'], mapDict(key) );
		catch e
			if strcmp(e.identifier,'MATLAB:Containers:Map:NoKey')
				warning(['Could not replace ',full,'. Leaving as is.']);
				val = '??';
			else
				rethrow e;
			end
		end
	end
end