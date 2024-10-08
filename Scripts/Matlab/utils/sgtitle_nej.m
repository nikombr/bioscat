function out = sgtitle(titleText, varargin)
    % Custom sgtitle function with LaTeX interpreter by default.
    % varargin allows for additional optional parameters if needed.

    % Check if 'Interpreter' is already specified
    if ~any(strcmpi(varargin, 'Interpreter'))
        varargin = [varargin, {'Interpreter', 'latex'}];
    end

    % Call the original sgtitle with modified arguments
    out = builtin('sgtitle', titleText, varargin{:});
end