%% INIT Reads the descriptor matrix in raw format and creates the internal data structure
function init(jac_val, jac_eqs, jac_vars, jac_struc)
global dif_eqs alg_eqs dif_states alg_states raw_vars fig S nbbus adf nbsync nbinj nbtwop bus_inj V eigenvals W eigenvals2 E Nx

%% Read values file
raw_val = importdata(jac_val) ;
% A Jacobian containing NaN or Inf is accepted silently by spconvert and turns
% the whole spectrum into NaN several steps later, so it is rejected here where
% the cause is still obvious. Octave returns a non-numeric value from
% importdata when the file holds NaN text, MATLAB returns the NaN itself.
if ~isnumeric(raw_val) || size(raw_val,2) < 3
    error('init:badValFile', ...
          ['%s is not a plain 3-column coordinate file. Non-numeric entries ' ...
           'such as NaN in the value column are the usual cause.'], jac_val);
end
if ~all(isfinite(raw_val(:,3)))
    error('init:nonFiniteJacobian', ...
          '%s holds %d non-finite entries (NaN or Inf) of %d; the spectrum would be all NaN.', ...
          jac_val, sum(~isfinite(raw_val(:,3))), size(raw_val,1));
end
S = spconvert(raw_val) ;
E = sparse(size(S,1),size(S,1)) ;

%% Read structure file if it exists
if ~isempty(jac_struc)
    fid = fopen (jac_struc);
    raw_struc = textscan(fid, '%d %d %d %d');
    fclose(fid);

    nbbus  = raw_struc{1}(1);
    nbsync = raw_struc{2}(1);
    nbinj  = raw_struc{3}(1);
    nbtwop = raw_struc{4}(1);

    adf = zeros(nbsync+nbinj+nbtwop+1,1);
    bus_inj = zeros(nbsync+nbinj+nbtwop+1,1);
    for idx = 2:numel(raw_struc{2})
        adf(idx-1) = raw_struc{2}(idx);
        bus_inj(idx-1) = raw_struc{3}(idx);
    end
end    

%% Read equations file
fid = fopen (jac_eqs);
raw_eqs = textscan(fid, '%d %s %s %s %s %d');
fclose(fid);

%% Read variables file
fid = fopen (jac_vars);
raw_vars = textscan(fid, '%d %s %s %s %s');
fclose(fid);

dif_states = []; % saves the places of differential states
alg_states = []; % saves the places of algebraic states

for idx = 1:numel(raw_vars{2})
    if char(raw_vars{2}(idx)) == 'd'
%         disp(strcat('element: ',num2str(idx),' is differential'))
        dif_states = [dif_states idx];
    else
        alg_states = [alg_states idx];
    end
end

dif_eqs=[]; % saves the places of differential equations
alg_eqs=[]; % saves the places of algebraic equations

gamma = sparse(numel(raw_eqs{6}),1);
for idx = 1:numel(raw_eqs{6})
    gamma(idx)=double(raw_eqs{6}(idx));
end

assignin('base', 'gamma', gamma);

% Put equations in order so that equation i is the derivative of state i
% This will make the matrix E of the descriptor (A,E) diagonal
idx = 1;
while idx <= numel(raw_eqs{6})
    if (gamma(idx) ~= 0) && (gamma(idx) ~= idx)
        % This equation should be in position gamma(idx) and not idx
        idx2 = gamma(idx);
        % disp(strcat('swapping: ',num2str(idx),' with ', num2str(idx2)))
        % We swap the matrix rows to put equation idx to the location idx2
        S([idx,idx2],:) = S([idx2,idx],:);
        % disp([num2str(idx) num2str(idx2)]);
        % We update gamma to reflect the swap
        gamma([idx,idx2]) = gamma([idx2,idx]);
        % We need to revisit new equation idx, since we don't know if this
        % is correct now. We only know that the one in the location idx2 is
        % now correct
        idx = idx - 1;
    end
    idx = idx +1;
end
assignin('base', 'gamma_rearranged', gamma);
assignin('base', 'S', S);

% Create the E matrix of the descriptor
for idx = 1:numel(raw_eqs{6})
    if gamma(idx) > 0
%         disp(strcat('element: ',num2str(idx),' is differential'))
        dif_eqs = [dif_eqs idx];
        E(idx,gamma(idx))=1;
    else
        alg_eqs = [alg_eqs idx];
    end
end

% Initialise the matrices.
V = [];
eigenvals = [];
eigenvals2 = [];
W = [];

Cdelta = sparse(length(dif_states),1);
Comega = sparse(length(dif_states),1);
B = sparse(length(dif_states),1);
for k = 1:length(dif_states)
    if strcmp(char(raw_vars{5}(dif_states(k))),'omega') 
        Cdelta(k)=1;
    elseif strcmp(char(raw_vars{5}(dif_states(k))),'delta')
        Comega(k)=1;
    elseif strcmp(char(raw_vars{5}(dif_states(k))),'vf')
        B(k)=1;
    end
end

Nx = numel(dif_eqs);
assignin('base', 'Cdelta', Cdelta);
assignin('base', 'Comega', Comega);
assignin('base', 'B', B);
assignin('base', 'E', E);
assignin('base', 'Nx', Nx);

close all
hold on
fig = 1;

fprintf('Number of diff-alg states=%d\n', numel(dif_states) + numel(alg_states));
fprintf('Number of diff states=%d\n', numel(dif_states));
end