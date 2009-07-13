classdef ProfileFunctions
    
	methods(Static)
 	
		function addedLine = addLine ( label, height, data, lineWidth, lineType, lineColor, lineCollection )
			%Adds a line to the created plot
			addedLine = plot(data, height, lineType, 'Color', lineColor, 'LineWidth', lineWidth);

			hold on
		end

		function addLegend ( lineCollection, textCollection )
			hold off

			legend( lineCollection, textCollection, 'Location', 'NorthEast' )
		end

		function setTitle ( graphTitle )
			hold off
		
			title(graphTitle)
		end

		function setAxisLabels ( xLabel, yLabel )
			hold off
		
			xlabel(xLabel)
			ylabel(yLabel)
		end
		
		function setAxis ( minVal, maxVal, startHeight, endHeight )
			hold off
		
			percentage = 95;
			
			% Normalize units to a factor of 1.
			percentage = percentage / 100;
			
			% The variable 'equiv_space' is the pure number of units of the variable 
			% being graphed to be displayed on either side of the vertical line if 
			% the variable is only equal to one value in the graph.  For example, if 
			% all values of thlm at every height for every data set are equal to 300 K,
			% then the plot will only show vertical lines at 300 K.  The number of
			% units (K, in this case) to the left or to the right (all white space)
			% will be determined by 'equiv_space'.  For this case, a value of 0.01
			% would make the graph range between 299.99 and 300.01.  The value of
			% 'equiv_space' must be a positive number.
			equiv_space = 0.01;
		
			% Extent of graph.
			xmin = min(minVal);
			xmax = max(maxVal);
			
			if ( xmax == xmin )
   				xmin = xmin - equiv_space;
   				xmax = xmax + equiv_space;
			else
   				xdiff = xmax - xmin;
   				xrange = xdiff / percentage;
   				xmedian = ( xmin + xmax ) / 2;
   				xmin = xmedian - xrange / 2;
   				xmax = xmedian + xrange / 2;
			end
		
			axis([ xmin xmax startHeight endHeight ])
		end	

    	end
    
end
