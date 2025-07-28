
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
    idx = event_obj.DataIndex;

    if isprop(fig, 'UserData') && iscell(fig.UserData)
        filenames = fig.UserData;
        if idx > 0 && idx <= numel(filenames)
            output_txt = {['File: ', filenames{idx}]};
        else
            output_txt = {'[Index out of bounds]'};
        end
    else
        output_txt = {'[No filenames found]'};
    end
end