
function output_txt = showFilename(event_obj, fig)
%SHOWFILENAMEF Custom data tip function to display file names.
%
%   OUTPUT_TXT = SHOWFILENAMEFROMFIGURE(EVENT_OBJ, FIG) returns a cell
%   array of text to be displayed in a data tip when the user clicks on a
%   point in a scatter plot. This function assumes that the full list of
%   file names corresponding to each data point is stored in the
%   'UserData' property of the figure FIG.
%
%   INPUTS:
%       event_obj - A handle to the data cursor event object (passed
%                   automatically by MATLAB's data tip system). Contains
%                   information about the clicked data point, including
%                   DataIndex and Position.
%
%       fig       - A handle to the figure object that contains the plot.
%                   The figure's UserData must contain a cell array of
%                   filenames (one per data point).
%
%   OUTPUT:
%       output_txt - A cell array of strings to be displayed in the data
%                    tip. In this case, it shows the filename associated
%                    with the clicked data point.
%
%   EXAMPLE USAGE:
%       fig = openfig('myplot.fig');
%       fig.UserData = filenames;
%       dcm = datacursormode(fig);
%       dcm.Enable = 'on';
%       dcm.UpdateFcn = @(src, event) showFilenameFromFigure(event, fig);
%
%   See also DATACURSORMODE, OPENFIG, USERDATA.

    % Find all scatter handles in plotting order
    scatterHandles = findobj(fig, 'Type', 'Scatter');
    scatterHandles = flipud(scatterHandles); % Ensure plotting order

    % Find which scatter object was clicked
    thisScatter = event_obj.Target;

    % Count how many points are before this object in plotting order
    offset = 0;
    for i = 1:length(scatterHandles)
        if scatterHandles(i) == thisScatter
            break;
        end
        offset = offset + numel(scatterHandles(i).XData);
    end

    % Global index = local index + offset
    globalIdx = event_obj.DataIndex + offset;

    filenames = fig.UserData;
    if globalIdx >= 1 && globalIdx <= numel(filenames)
        output_txt = {['File: ', filenames{globalIdx}]};
    else
        output_txt = {'[Index out of bounds]'};
    end
end
