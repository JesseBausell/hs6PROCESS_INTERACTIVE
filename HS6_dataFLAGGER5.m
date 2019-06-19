function [flaggedROWS, HS6_data2] = HS6_dataFLAGGER5(HS6_data,lambdA,flaggedROWS)
% HS6_dataFLAGGER5
% Jesse Bausell
% November 18, 2016
%
% HS6_dataFLAGGER5 allows user to determine interactively which data to flag
% for later removal and which data to keep. User's selections will be saved
% to be used down the road in the actual deselecting program. The program
% also lets the user discard sections of data that are inaccurate, by
% specifying a range of depths to eliminate, but only at the beginning of
% the program. Upgraded from HS6_dataFLAGGER3. This version automatically
% flags HS6 values that are outside of the user-prescribed range.
%
%
% Inputs:
% HS6_data - sigma-corrected HS6 data to be flagged
% LambdA - the wavelength of the HS6 data
% flaggedROWS - the rows already flagged (this is used 
%
% Outputs:
% flaggedROWS - row indices that have been flagged by the user.
% HS6_data - HS6 data NOT DISCARDED by the user
%% 1. Prepare for the plotting

% Determine the number of HS6 Channels and make plotting variables that fit
% them.

        close all; clc;
        % Closes all plots and clears the command window.
        
        HS6_data2 = HS6_data;
        % Create a dublicate HS6 data set.
        
        [l,w] = size(HS6_data);
        % Create length and width variables for HS6_data.
        
        deptH = 1:l;
        % Create depth index array to be used in subsequent analysis.
        
        if isequal(length(lambdA),6)
            % This if statement determines whether we are dealing with 6 or
            % 8 bb channels in the HS6. This line deals with 6.

            coloR = [.6 0 1;0 0 1;0 .6 0;.8 .8 0;.8 0 0;0 0 0];
            % Color codes for six HS6 channels
            
            legendTEXT = {[num2str(lambdA(1)) 'nm '], [num2str(lambdA(2)) 'nm'], ...
                [num2str(lambdA(3)) 'nm'], [num2str(lambdA(4)) 'nm'], ...
                [num2str(lambdA(5)) 'nm'], [num2str(lambdA(6)) 'nm'], 'flagged'};
            % Legend names (incorporating lambdas) for six HS6 channels

        else
            % IF there are 8 bb channels in the HS6, use these options instead 

            coloR = [.6 0 1; 0 0 1; 0 .8 .8; 0 .6 0; .8 .8 0; .8 0 0; 1 .5 .3; 0 0 0];
            % Color codes for eight HS6 channels

            legendTEXT = {[num2str(lambdA(1)) 'nm '], [num2str(lambdA(2)) 'nm'], ...
                [num2str(lambdA(3)) 'nm'], [num2str(lambdA(4)) 'nm'], ...
                [num2str(lambdA(5)) 'nm'], [num2str(lambdA(6)) 'nm'], ...
                [num2str(lambdA(7)) 'nm'], [num2str(lambdA(8)) 'nm'],'flagged'};
            % Legend names (incorporating lambdas) for eight HS6 channels

        end

        
        if isnan(flaggedROWS)
            % If flaggedROWS input is NaN...  
            flaggedROWS = []; % Change this variable to empty
            locK = 0; % This variable will allow user to permanently discard data.    
        else             
            locK = 1; % This variable will allow user to permanently discard data.
        end
    
        flagged_IND = NaN; 
        % creates a variable in which to put newly flagged HS6 rows
        
         levER = 0; 
         % This variable determines which plotting mode is used:
         % 0 = displays all depth profiles & flagged points
         % 1 = displays depth profiles with newly-selected flagged point

%% 2. Next we do the dirty work by deciding which data to discard, flag and keep.

    while 1
        
        %figure('Color',[.8 .8 .8]);
        % sets a figure with a grey background color.
        
        % Creates an ongoing loop with which to flag questionable data.
         
%          flagged_IND = NaN; 
         % creates a variable in which to put newly flagged HS6 rows

         depth_boTTom = -.2*max(deptH);
         depth_tOp = max(deptH)+2;
         HS6_toP = max(HS6_data(:))*1.1;
         HS6_boTTom = min(HS6_data(:))-.001;
         % These variables will determine the range of the subsequent HS6 depth
         % profiles.

        hold on; %Plot multiple things on the same graph.
        
        for ii = 1:length(lambdA)
            %Plot all HS6 channels one at a time. 
            
            handLE(ii).HS6 = plot(deptH,HS6_data(:,ii),'Color',coloR(ii,:));
   
        end
    
         if isequal(levER,0)
                    % This if statement differentiates between two paths,
                    % responsible for two different graphs. Path 1 (this
                    % path), will plot HS6 depth profiles along with
                    % previously-flagged data. levER = 0.

            if ~isempty(flaggedROWS)
                % If the user has previously flagged questionable data
                
                for jj = 1:length(flaggedROWS)
                % Plot flagged data on top of profiles (as black stars) one
                % row at a time.
                
                    h = scatter(deptH(flaggedROWS(jj)*ones(1,length(lambdA))),HS6_data(flaggedROWS(jj),:),'w*');
                    % Plot points on top of depth profiles

                end    
            
            else
                % If this is the first time that the user is flagging
                % questionable data.
                
                    h = scatter(NaN,NaN,'w*');
                    % This is a dummy plot. It doesn't actually change the
                    % appearance, but it allows us to use the same legend.
                
            end

                plot([depth_boTTom depth_tOp],[0 0],'k'); plot([0 0],[HS6_boTTom HS6_toP],'k');
                % This line of code plots reference lines to show a depth
                % of zero m and a HS6 of zero m^-1. The purpose of these
                % are to give users a reference with which to flag bad
                % data.
                
                axis([depth_boTTom depth_tOp HS6_boTTom HS6_toP]);
                view(-90,90);
                set(gca,'XAxisLocation','bottom','YAxisLocation','right', ...
                    'ydir','reverse','xdir','reverse','FontSize',10,'XTick', ...
                    [0:100:ceil(max(deptH))],'YTick', ...
                    [min(HS6_data(:)):abs(max(HS6_data(:))-min(HS6_data(:)))/5:max(HS6_data(:))],'Color',[.8 .8 .8]);
                xlabel('Depth (m)','FontSize',16);
                ylabel('Backscatter','FontSize',16);
                % These code lines flip the profile to make it vertical,
                % and 'decorate' the figure properly.
            
                if isequal(length(lambdA),6)
                    % This if statement provides the proper legend: 6
                    % channels vs. 8 channels
                    
                    LEG = legend([handLE(1).HS6 handLE(2).HS6 handLE(3).HS6 ...
                        handLE(4).HS6 handLE(5).HS6 handLE(6).HS6 h(1)]...
                        ,legendTEXT,'Location','northeast');
                    % Legend for 6 channels
                    
                    set(LEG,'FontSize',12) % set font to 12

                else
                    
                     LEG = legend([handLE(1).HS6 handLE(2).HS6 handLE(3).HS6 ...
                         handLE(4).HS6 handLE(5).HS6 handLE(6).HS6 handLE(7).HS6...
                         handLE(8).HS6 h(1)],legendTEXT,'Location','northeast');
                    % Legend for 8 channels
                    
                     set(LEG,'FontSize',12) % set font to 12.

                end
                
                if isequal(locK,0)
                    % This if statement offers two different paths
                    % depending on whether or not the user wants to discard
                    % data. locK = 0 indicates that the user a) has not yet
                    % decided if he/she wants to discard data or b) has
                    % indicated an interest to discard data. Once user
                    % decides not to discard data, this path is closed
                    % permanently.
                    
                        while 1
                            % This is a while loop within a while loop. It allows
                            % the user to limit the boundaries of the HS6 data if
                            % he/she feels that they are very far off. 

                            keY5 = input('Create HS6 limit? (y/n)','s');
                            % Asks user if he/she wants to limit data.

                            if strcmpi('y',keY5)
                                %If the user says 'yes', he/she will be prompted to
                                %choose an upper and lower boundary.
                                hold on;
                                uppER = input('Select upper HS6 limit: '); %upper boundary selection
                                lowER = input('Select lower HS6 limit: '); %lower boundary selection

                                plot([depth_boTTom depth_tOp],[uppER uppER],'k','LineWidth',1.5);
                                plot([depth_boTTom depth_tOp],[lowER lowER],'k','LineWidth',1.5);
                                % Plots upper and lower boundary to give user a
                                % visual on which data will be axed.

                                keY6 = input('Accept these limits? (y/n)','s');
                                % Lets user sign off on his/her choice

                                    if strcmpi(keY6,'y')
                                        % If user confirms his/her selection

                                        limit_IND = (lowER > HS6_data | HS6_data > uppER);
                                        % Finds all the falling outside the user's
                                        % new range.

                                        HS6_data(limit_IND) = NaN;
                                        % Replaces all ouside data with NaNs
                                        
                                        [L_nan,w_nan] = find(isnan(HS6_data));
                                        % Find all nans in the HS6 data (nans denote data that was
                                        % outside of our acceptable range.

                                        flaggedROWS = L_nan;
                                        % Add nan'd values to the flaggedROWS.

                                        flaggedROWS = sort(flaggedROWS);
                                        % Puts row indices in order

                                        indexROWS = ones(length(flaggedROWS),1);
                                        % Gives us a dummy array with which to match up to
                                        % flaggedROWS variable.

                                        for ii = 2:length(flaggedROWS)
                                            % Go through all flaggedROWS

                                            if isequal(flaggedROWS(ii-1),flaggedROWS(ii))
                                                % If a row number repeats itself, mark it as NaN.

                                                indexROWS(ii) = NaN;
                                            end

                                        end

                                        flaggedROWS = flaggedROWS(~isnan(indexROWS));
                                            % Extract the nans from the flaggedROWS, leaving only
                                            % the actual values.

                                        
                                        close all; %closes the graph
                                        break %breaks the inside while loop

                                    end

                            elseif strcmpi('n',keY5)
                                        % If the user decides he/she doesn't want
                                        % to select any more data limits.

                                        locK = 1; 
                                        %changes the lock variable to 1 (which
                                        %opens up the next part of the program).
                                        close all; %closes the plot
                                        break %breaks the inner while loop

                            else

                                clc; %clears the command window
                                disp('Invalid entry. Must select "y" or "n"');
                                % Informs user that he/she selected an invalid
                                % entry.

                            end
                            
                        end      
                        
                elseif isequal(locK,1)
                    % This option, locK = 1, indicates that the user is
                    % finished discarding data or has expressed a desire
                    % NOT to discard data.We now move on to
                    % flagging individual data points.
                
                    keY = input('Flag additional HS6 readings? (y/n): ','s');
                    % Lets user decide whether to flag a data point that
                    % may be questionable.

                    if strcmpi('n',keY)
                        % No additional data points are flagged. Function
                        % is halted.

                        clc; close all; break
                        % clear command window, close figures, and break
                        % the while loop

                    elseif strcmpi('y',keY)
                        % User decides he/she wants to flag a data point. 

                        levER = 1;
                        % Lever changes to 1. This will open up a new
                        % figure where user can make sure he/she selected
                        % the correct data point. 

                        deptH_KEY = input('Enter depth of data you want to flag: ');
                        % User enters the depth of the data row he/she
                        % wishes to flag.

                        deptH_DIFF = abs(deptH-deptH_KEY);
                        % Finds the actual depth that is closest to the
                        % value he/she entered.

                        flagged_IND = find(deptH_DIFF == min(deptH_DIFF));
                        flagged_IND = flagged_IND(1);
                        % Generates the index of the depth of interest

                        close all; %closes figure. While loop begins a new.
                    end
            
                end

        elseif isequal(levER,1)
                        % This path signifies that the user has selected a
                        % new HS6 row to flag (potentially). This path will
                        % graph the HS6 depth profiles and put a black star
                        % on the flagged data point to let the user decide
                        % whether to permanently flag it.
            
            h = scatter(deptH(flagged_IND)*ones(1,length(lambdA)),HS6_data(flagged_IND,:),'w*');
            % This line of code will plot the row of questionable data as a
            % one row scatter plot.
                
                axis([depth_boTTom depth_tOp HS6_boTTom HS6_toP]);
                view(-90,90);
                set(gca,'XAxisLocation','bottom','YAxisLocation','right', ...
                    'ydir','reverse','xdir','reverse','FontSize',10,'XTick', ...
                    [0:100:ceil(max(deptH))],'YTick', ...
                    [min(HS6_data(:)):abs(max(HS6_data(:))-min(HS6_data(:)))/5:max(HS6_data(:))], ...
                    'Color',[.8 .8 .8]);
                xlabel('Depth (m)','FontSize',16);
                ylabel('Backscatter','FontSize',16);
                % These code lines flip the profile to make it vertical,
                % and 'decorate' the figure properly.
                
                if isequal(length(lambdA),6)
                    % This if statement provides the proper legend: 6
                    % channels vs. 8 channels
                    
                    LEG = legend([handLE(1).HS6 handLE(2).HS6 handLE(3).HS6 ...
                        handLE(4).HS6 handLE(5).HS6 handLE(6).HS6 h(1)]...
                        ,legendTEXT,'Location','northeast');
                    % Legend for 6 channels
                    
                    set(LEG,'FontSize',12) % set font to 12

                else
                    
                     LEG = legend([handLE(1).HS6 handLE(2).HS6 handLE(3).HS6 ...
                         handLE(4).HS6 handLE(5).HS6 handLE(6).HS6 handLE(7).HS6...
                         handLE(8).HS6 h(1)],legendTEXT,'Location','northeast');
                    % Legend for 8 channels
                    
                     set(LEG,'FontSize',12) % set font to 12.

                end

         
                while 1
                    % This while loop within a while loop sets up the
                    % foundation for the final question: whether to
                    % actually flag the data.
                    
                    keY2 = input('Flag this HS6 row? (y/n) ','s');
                    % Asks user if he/she is sure they want to flag the
                    % data.

                    if strcmpi('y',keY2)
                        % If user decides that they want to flag the data.
                            
                            flaggedROWS = [flaggedROWS; flagged_IND];
                            % Adds the flagged row to the permanent array
                            % of flagged rows.
                            
                            flaggedROWS = sort(flaggedROWS);
                            % sorts flagged rows in proper order
                            % (ascending)
                            
                            flagged_IND = NaN; levER = 0;
                            % Resets key variables back to their original
                            % values. flagged_IND carries newly flagged
                            % data and levER puts program back on the
                            % original track, to plot depth profiles with
                            % ALL flagged points next time arround.
                            
                            break %breaks the nested while loop.
                            
                    elseif strcmpi('n',keY2)
                        % If user decides NOT to flag the data.

                        flagged_IND = NaN; levER = 0;
                        % Key variables reset to original values (see above
                        % for explanation)
                            
                            break % breaks the nested while loop.
                        
                    end
                
                end

          end
        
    end
    
    
%% 3. Finally we decide whether to keep modified HS6 data.

    while 1 
        % Final while loop to decide whether to keep or discard 'edited'
        % HS6 data.
        
        keY7 = input('Accept modified HS6 data? (y/n)','s');
        % Question prompting user on whether to keep or discard HS6 data.

            if strcmpi(keY7,'y')
                % If user answers "yes", the modified data will become the
                % primary HS6 data set to be used in
                % HydroScat_Processor_DOXARAN.
                
                HS6_data2 = HS6_data; % change data set
                break %break the while loop/end the function

            elseif strcmpi(keY7,'n')
                 
                [L_nan,w_nan] = find(isnan(HS6_data));
                % Find all nans in the HS6 data (nans denote data that was
                % outside of our acceptable range.
                
                flaggedROWS = sort(flaggedROWS);
                % Puts row indices in order
                
                indexROWS = ones(length(flaggedROWS),1);
                % Gives us a dummy array with which to match up to
                % flaggedROWS variable.
                
                for ii = 2:length(flaggedROWS)
                    % Go through all flaggedROWS
                    
                    if isequal(flaggedROWS(ii-1),flaggedROWS(ii))
                        % If a row number repeats itself, mark it as NaN.
                        
                        indexROWS(ii) = NaN;
                    end
                    
                end
                
                flaggedROWS = flaggedROWS(~isnan(indexROWS));
                    % Extract the nans from the flaggedROWS, leaving only
                    % the actual values.
                
                break % break the while loop/end the function
                
            else
                
                disp('Invalid entry. Try again');
                % Inform user that he/she's a dumbass.

            end
    
    end
    
end
