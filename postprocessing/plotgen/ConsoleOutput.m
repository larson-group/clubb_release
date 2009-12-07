classdef ConsoleOutput
    
	methods(Static)
 	
		function warning ( text )
		    unix(['./console_output.pl ' '-w ' '"' text '"']);
        end

        function message ( text )
            unix(['./console_output.pl ' '"' text '"']);
        end

        function severe ( text )
            unix(['./console_output.pl ' '-s ' '"' text '"']);
        end
    end
end
