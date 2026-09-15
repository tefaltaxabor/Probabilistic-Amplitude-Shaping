function no_errors_legend(ax, label)
% NO_ERRORS_LEGEND  One legend entry explaining the src.mark_no_errors markers.
%
%   src.no_errors_legend(ax)
%   src.no_errors_legend(ax, label)
%
%   Adds an invisible dummy series (grey open downward triangle) whose only
%   purpose is the legend entry.

    if nargin < 2
        label = 'no errors (upper bound 1/N)';
    end
    plot(ax, NaN, NaN, 'v', 'Color', [0.35 0.35 0.35], 'MarkerSize', 7, ...
         'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', label);
end
