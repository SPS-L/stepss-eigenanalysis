%% CAPTURE_GOLDEN Writes reference outputs of the QZ path as plain text
%
%  Copyright 2026 Sustainable Power Systems Lab (SPSL)
%
%  Licensed under the Apache License, Version 2.0 (the 'License');
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0
%
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an 'AS IS' BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%
%  Usage:
%  capture_golden('jac', 'golden')
%  capture_golden('jac', 'golden', 'jac_struc.dat')
%
%  Reads <prefix>_val.dat, <prefix>_eqs.dat and <prefix>_var.dat, runs the
%  reduction and the dense eigensolve, and writes four text files into outdir.
%  This is the reference the Fortran and Python implementations are validated
%  against, so it deliberately avoids analyze_results and loop_analysis: no
%  plotting, no interaction, no base-workspace saving.
%
%  Runs identically under MATLAB and Octave 8. From a shell:
%    octave --no-window-system --quiet --eval \
%      "addpath('..','../scripts'); capture_golden('test','../golden')"
%
%  Outputs, all 1-based and full double precision:
%    <name>_Asys.txt    i j value, one nonzero of the state matrix per line
%    <name>_eigs.txt    real imag, one eigenvalue per line
%    <name>_pf.txt      dense participation factors, rows = differential
%                       states, columns = modes, normalised so each column
%                       peaks at 1
%    <name>_states.txt  index family device variable, one differential state
%                       per line, giving the row labels of _pf.txt
%
%  Modes are sorted by descending real part, then descending imaginary part,
%  so the files are reproducible. Note that eigenvectors are NOT captured:
%  their scale and phase are arbitrary, and for degenerate eigenvalues they
%  are not even unique, so they cannot serve as a reference. Participation
%  factors are only meaningful for simple eigenvalues for the same reason.

function capture_golden(prefix, outdir, jac_struc)

global dif_states raw_vars Jdyn Jdyn_sp

if nargin < 3 || isempty(jac_struc)
    jac_struc = '';
end
if nargin < 2 || isempty(outdir)
    error('capture_golden requires a prefix and an output directory.');
end

[~, name] = fileparts(prefix);

% init rejects a non-finite Jacobian, so a NaN export fails here with a clear
% message rather than producing an all-NaN spectrum.
init([prefix '_val.dat'], [prefix '_eqs.dat'], [prefix '_var.dat'], jac_struc);
calc_Jdyn();

[V, D, W] = eig(Jdyn);
lam = diag(D);

%% Deterministic order: descending real part, then descending imaginary part
[~, ord] = sortrows([real(lam), imag(lam)], [-1 -2]);
lam = lam(ord); V = V(:,ord); W = W(:,ord);

%% Participation factors, as analyze_results.m computes them
P = abs(W .* V);
for i = 1:size(P,2)
    P(:,i) = P(:,i) / max(P(:,i));
end

%% Residual check: independent of any reference, so it must hold here
res = norm(Jdyn*V - V*diag(lam), 'fro') / norm(Jdyn, 'fro');
fprintf('%s: Nx=%d, relative eigenpair residual %.3e\n', name, numel(lam), res);
if res > 1e-10
    error('Eigenpair residual %.3e is too large to use as a reference.', res);
end

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

[ri, ci, vi] = find(Jdyn_sp);
write_rows([outdir filesep name '_Asys.txt'], '%d %d %.17e\n', [ri(:), ci(:), vi(:)]);
write_rows([outdir filesep name '_eigs.txt'], '%.17e %.17e\n', [real(lam), imag(lam)]);

fid = fopen([outdir filesep name '_pf.txt'], 'w');
fprintf(fid, [repmat('%.17e ', 1, size(P,2)) '\n'], P.');
fclose(fid);

fid = fopen([outdir filesep name '_states.txt'], 'w');
for k = 1:numel(dif_states)
    fprintf(fid, '%d %s %s %s\n', dif_states(k), ...
            char(raw_vars{3}(dif_states(k))), ...
            char(raw_vars{4}(dif_states(k))), ...
            char(raw_vars{5}(dif_states(k))));
end
fclose(fid);

fprintf('%s: wrote 4 golden files to %s\n', name, outdir);
end

function write_rows(fname, fmt, data)
fid = fopen(fname, 'w');
fprintf(fid, fmt, data.');
fclose(fid);
end
