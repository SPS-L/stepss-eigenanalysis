%% SSA Receives the descriptor form of an eigenvalue problem and calculates the eigenvalues/eigenvectors with different methods
%  
%  Copyright 2024 Sustainable Power Systems Lab (SPSL)
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
%  ssa('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat',real_limit,damp_ratio,method)
%  ssa('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat',real_limit,damp_ratio,[])
%  ssa('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat',real_limit,[],[])
%  ssa('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat',[],[],[])
%
%  jac_val : has the values of the matrix in coordinate format
%  jac_eqs : has the description of the equations, mainly if they are
%  differential or algebraic
%  jac_val : has the description of the variables, mainly if they are
%  differential or algebraic
%  jac_struc : decomposed structure of power system
%  real_limit : the real number above which the eigenvalue is considered as
%  dominant (optional, default=-inf)
%  damp_ratio : Dumping ratio above which the eigenvalue is considered as
%  dominant (optional, default=1.0)
%  Runs under MATLAB R2016a or later and under GNU Octave 8. Only the 'QZ'
%  method is supported on Octave: 'ARP' calls eigs() for every finite
%  eigenvalue of a pencil with singular E, which Octave's ARPACK refuses with
%  "Could not build an Arnoldi factorization".
%
%  method: The default method is 'QZ' through eig() which requires to make an
%  ellimination of the algebraic variables. If you pass 'ARP' then it will
%  use the descriptor method with the sparse matrices (A,E) and the
%  libraries provided by ARP (Krylov-Schur Algorithm) through eigs(). (optional, default='QZ')

function ssa(jac_val,jac_eqs,jac_vars,jac_struc,real_limit,damp_ratio,method)

global sigma analysis numEig

%% Check and correct arguments
if nargin < 3
    error('You need to give at least 3 arguments. Write ''help ssa''.');
elseif nargin > 7
    error('ssa requires at most 7 arguments. Write ''help ssa''.');
end


if ~exist('method','var') || isempty(method)
    analysis = 'QZ' ; 
else 
    analysis = method ; 
end
if ~(strcmp(analysis,'QZ') || strcmp(analysis,'ARP'))
    error('Method should be ''QZ'' or ''ARP''.');
end

if ~exist('jac_struc','var') || isempty(jac_struc)
  jac_struc='';
end

if ~exist('real_limit','var') || isempty(real_limit)
  real_limit = -inf;
end

if ~exist('damp_ratio','var') || isempty(damp_ratio)
  damp_ratio = 1.0;
end


%% Initialize data structure
init_tmr = tic;
init(jac_val, jac_eqs, jac_vars, jac_struc);
fprintf('Initialization done in %.3f seconds.\n\n', toc(init_tmr));

%% Calculate eigenvalues and eigenvectors using eigs (ARPACK Arnoldi method)
% Eigs_tmr=tic;
% analysis = 'IRA';
% numEig = 10 ;
% real_part = -0.1;
% while real_part > -1.0
%     freq = sqrt((real_part^2-damp_ratio^2*real_part^2)/damp_ratio^2)/(2.0*pi);
%     while freq < 2.0
%         sigma = real_part+2.0*pi*freq*1i;
% %         sigma = 0.0 + 6.28*1i;
% %         freq=6.28/(2.0*pi);
%         fprintf('Sigma= %.2f + i %.2f (freq= %.2f).\n',real(sigma),imag(sigma),freq);
%         eigenvals_eigs();
%         sigma = real_part-2.0*pi*freq*1i;
% %         sigma = -0.5 - 2.0*pi*0.5*1i;
%         fprintf('Sigma= %.2f + i %.2f (freq= %.2f).\n',real(sigma),imag(sigma),freq);
%         eigenvals_eigs();
%         freq = freq +0.1;
%     end
%     real_part = real_part - 0.1;
% end
% sigma = 0.1*1i;
% eigenvals_eigs();
% analyze_results(real_limit,damp_ratio);
% fprintf('Total time spent in eigs %.3f seconds.\n\n',toc(Eigs_tmr));

%% Calculate eigenvalues and eigenvectors using JDQR
% Eigs_tmr=tic;
% eigenvals_jdqr();
% fprintf('Total time spent in eigs %.3f seconds.\n\n',toc(Eigs_tmr));

%% Calculate eigenvalues and eigenvectors
if strcmp(analysis,'QZ')
    eigenvals_eig();
elseif strcmp(analysis,'ARP')
    eigenvals_eig_descr();
end
analyze_results(real_limit, damp_ratio);
% '-v7' is given explicitly so MATLAB and Octave both write a MAT file of the
% same name and format. Octave's bare save() would write its own format to a
% file with no extension.
evalin('base', 'save(''modal_reduction.mat'', ''-v7'')');
evalin('base', 'savefig(''eigs'')');
loop_analysis();
% fprintf('\nExecute ssa to rerun the analysis.\n');
end
