function save_result(fig, name, outdir, data)
% SAVE_RESULT  Persist one showcase panel: .fig + .pdf (+ .mat with the data).
%
%   src.save_result(fig, name, outdir, data)
%
%   Follows the project convention of writing both the MATLAB figure and the
%   raw data into results/, so every plot can be regenerated or restyled for
%   the thesis without re-running the Monte Carlo. Safe to call headless
%   (matlab -batch on the cluster): figures are created with 'Visible','off'
%   by the callers and exportgraphics does not need a display.
%
%   fig    : figure handle (or [] to save data only)
%   name   : base filename, without extension
%   outdir : destination folder (default 'results')
%   data   : optional struct; its fields are saved as variables in name.mat

    if nargin < 3 || isempty(outdir), outdir = 'results'; end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    base = fullfile(outdir, name);

    if ~isempty(fig) && isgraphics(fig)
        savefig(fig, [base '.fig']);
        try
            exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
        catch err
            warning('save_result:pdf', 'PDF export failed for %s: %s', name, err.message);
        end
    end

    if nargin >= 4 && ~isempty(data)
        save([base '.mat'], '-struct', 'data');
    end

    fprintf('  saved -> %s.{fig,pdf,mat}\n', base);
end
