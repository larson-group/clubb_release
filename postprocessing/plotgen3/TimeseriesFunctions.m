classdef TimeseriesFunctions
    
	methods(Static)
 	
		function addedLine = addLine ( label, time, data, lineWidth, lineType, lineColor, lineCollection )
			%Adds a line to the created plot
			addedLine = plot(time, data, lineType, 'Color', lineColor, 'LineWidth', lineWidth);

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
		
		function setAxis ( minVal, maxVal, startTime, endTime )
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
			ymin = min(minVal);
			ymax = max(maxVal);
			
			if ( ymax == ymin )
   				ymin = ymin - equiv_space;
   				ymax = ymax + equiv_space;
			else
   				ydiff = ymax - ymin;
   				yrange = ydiff / percentage;
   				ymedian = ( ymin + ymax ) / 2;
   				ymin = ymedian - yrange / 2;
   				ymax = ymedian + yrange / 2;
			end
		
			axis([ startTime endTime ymin ymax ])
		end	

    	end
    
end
