function any0 = mark_no_errors(ax, x, y, yRes, color)
% MARK_NO_ERRORS  Draw zero-error Monte-Carlo points as upper bounds.
%
%   any0 = src.mark_no_errors(ax, x, y, yRes, color)
%
%   A simulated point with no errors does not have error rate 0 (which a log
%   axis cannot show and nan_zeros() hides), it has error rate < yRes, the
%   resolution of that point (1 / codewords or bits simulated). Those points
%   are drawn as open downward triangles at yRes in the colour of the curve,
%   not joined to it. Returns true if any such point was drawn, so the caller
%   can add a single legend entry with src.no_errors_legend.

    z = (y(:).' == 0);
    any0 = any(z);
    if ~any0, return; end
    x = x(:).';
    plot(ax, x(z), repmat(yRes, 1, nnz(z)), 'v', 'Color', color, ...
         'MarkerSize', 7, 'LineWidth', 1.2, 'LineStyle', 'none', ...
         'HandleVisibility', 'off');
end
