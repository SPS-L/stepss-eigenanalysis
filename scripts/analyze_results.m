%% ANALYZE_RESULTS Analyzes the results of the eigenvalue/eigenvector computation
function analyze_results(real_limit,damp_ratio)
global dif_states P dom_V fig dom_eigenvals DF eigenvals dom_W V W eigenvals2 analysis

% assignin('base', strcat('eigenvals', analysis), eigenvals);
% assignin('base', strcat('eigenvals2', analysis), eigenvals2);
% assignin('base', strcat('V', analysis), V);
% assignin('base', strcat('W', analysis), W);
% return

if strcmp(analysis,'IRA')
    
    tol = 10^-8;
    
    % remove duplicate values
    [sorted_eigenvals, sorted_V] = cplx_unique(eigenvals,V,tol);
    [sorted_eigenvals2, sorted_W] = cplx_unique(eigenvals2,W,tol);
    disp(sorted_eigenvals);
    disp(sorted_eigenvals2);
    
    % match right to left eigenvectors through the eigenvalues
    m = length(sorted_eigenvals);
    n = length(sorted_eigenvals2);
    kronp = transpose(reshape(abs(imag(kron(sorted_eigenvals, conj(sorted_eigenvals2))))<tol,n,m));
    disp(kronp);
    idxA = zeros(m,1);
    for i = 1:m
        for j=1:n
            if (kronp(i,j))
                idxA(i) = i; % it exists in the other list
            end
        end
    end
    
    idxA = idxA(idxA~=0);
    sorted_eigenvals = sorted_eigenvals(idxA); % remove eigenvalues not the same in both lists
    sorted_V=sorted_V(:, idxA); % remove eigenvectors not the same
    disp(sorted_eigenvals);
    
    m = length(sorted_eigenvals);
    n = length(sorted_eigenvals2);
    kronp = transpose(reshape(abs(imag(kron(sorted_eigenvals, conj(sorted_eigenvals2))))<tol,n,m));
    idxB = zeros(n,1);
    for i = 1:m
        for j=1:n
            if (kronp(i,j))
                idxA(i) = i;
                idxB(j) = i; % correspondence to other list
            end
        end
    end
    
    idxB = idxB(idxB~=0); % remove the ones that do not exist
    sorted_eigenvals2=sorted_eigenvals2(idxB);
    sorted_W=sorted_W(:,idxB);
    
    m = length(sorted_eigenvals);
    n = length(sorted_eigenvals2);
    if (m ~= n)
        error('At this point the left and right eigenvectors should be matching')
    end

elseif strcmp(analysis,'ARP')
    
    tol = 10^-8;
    
    % remove duplicate values
    [sorted_eigenvals, sorted_V] = cplx_unique(eigenvals,V,tol);
    [sorted_eigenvals2, sorted_W] = cplx_unique(eigenvals2,W,tol);
    
else
    [foo, sortidx] = sort(real(eigenvals),'descend');
    sorted_eigenvals = eigenvals(sortidx);
    sorted_V = V(:,sortidx);
    sorted_W = W(:,sortidx);
end

%% Compute dominant eigenvalues and dumping factors

dom_eig_num = zeros(length(sorted_eigenvals)); % addresses of dominant eigenvalues
DF = zeros(length(sorted_eigenvals)); % Dumping Factors
for idx = 1:length(sorted_eigenvals)
    a=real(sorted_eigenvals(idx));
    b=imag(sorted_eigenvals(idx));
    DF(idx)=-a/sqrt(a^2+b^2);
    if a > real_limit % && b>=0.0 % save only dominant eigenvalues and only the upper eigenvalue in the complex plane
        dom_eig_num(idx)=idx;
    end
end

dom_eig_num = dom_eig_num(dom_eig_num~=0);

if ~any(dom_eig_num)
    disp('No dominant eigenvalues with the criteria asked.')
    return
end

dom_eigenvals=sorted_eigenvals(dom_eig_num); % saves the dominant eigenvalues
if strcmp(analysis,'IRA')
    sorted_eigenvals2=sorted_eigenvals2(dom_eig_num); % saves the dominant eigenvalues
end
DF = DF(dom_eig_num);
dom_V = sorted_V(:,dom_eig_num); % saves the respective right eigenvector of dominant eigenvalues
dom_W = sorted_W(:,dom_eig_num); % saves the respective left eigenvector of dominant eigenvalues

assignin('base', strcat('eigenvals', analysis), dom_eigenvals);
if strcmp(analysis,'IRA')
    assignin('base', strcat('eigenvals2', analysis), sorted_eigenvals2);
end
assignin('base', strcat('V_', analysis), dom_V);
assignin('base', strcat('W_', analysis), dom_W);

%% Compute and normalize participation factors
P=zeros(length(dif_states),length(dom_eigenvals)); % Participation Factor matrix
for i = 1:length(dom_eigenvals)
    for k = 1:length(dif_states)
        P(k,i)=abs(dom_W(k,i)*dom_V(k,i)); % compute participation factor
    end
end

%  Normalize Participation Factor matrix
for i = 1:length(dom_eigenvals)
    pmax = max(P(:,i));
    %     disp(raw_vars{5}(dif_states(idx)))
    for k = 1:length(dif_states)
        P(k,i)=P(k,i)/pmax;
    end
end

%% Display Dominant eigenvalues
disp(strcat('Dominant eigenvalues: ', analysis))
for i = 1:length(dom_eigenvals)
    fprintf('Eigenvalue %3d:\t%12.6f+%12.6fi\tDF: %f\tf= %f\n',i,real(dom_eigenvals(i)),imag(dom_eigenvals(i)),DF(i),abs(imag(dom_eigenvals(i)))/(2*pi));
end
% assignin('base', 'qz_nordic', dom_eigenvals);

%% Plot dominant eigenvalues
figure(fig)
if strcmp(analysis,'IRA')
    scatter(real(dom_eigenvals),imag(dom_eigenvals),'filled', 'd', 'DisplayName', analysis);
else
    scatter(real(dom_eigenvals),imag(dom_eigenvals), 'DisplayName', analysis);
end
legend('-DynamicLegend');

% if strcmp(analysis,'IRA')
%     plot(real(dom_eigenvals),imag(dom_eigenvals),'r*')
% else
%     plot(real(dom_eigenvals),imag(dom_eigenvals),'bd')
%     legend('Arnoldi method','QZ method')
% end
% legend(strcat(analysis, ' method'))

xlabel('Real')
ylabel('Imaginary')
t1 = 'Dominant Eigenvalues';
title(t1)
xmin=min(min(real(dom_eigenvals)))-1;
xmax=max(max(real(dom_eigenvals)));
if xmax<0
    xmax=0 ;
end
ymin=-0.1 ;
ymax=max(max(imag(dom_eigenvals)))+1;
 
xlim([xmin xmax]);
ylim([ymin ymax]);
 
lh = line([real_limit,real_limit],[ymin,ymax],'Color','g');
set(lh, 'DisplayName', 'Dominant Eigenvalues Real Limit');
lh = line([0,0],[ymin,ymax],'Color','r');
set(lh, 'DisplayName', 'Imaginary Axis');

if ymin < xmin
    lh = line([0,ymin],[0,ymin*tan(acos(damp_ratio))],'Color','b');
    set(lh, 'DisplayName', 'DF Limit');
    lh = line([0,ymin],[0,ymax*tan(acos(damp_ratio))],'Color','b');
    set(lh, 'DisplayName', 'DF Limit');
else
    lh = line([0,xmin],[0,xmin*tan(acos(damp_ratio))],'Color','b');
    set(lh, 'DisplayName', 'DF Limit');
    lh = line([0,xmin],[0,-xmin*tan(acos(damp_ratio))],'Color','b');
    set(lh, 'DisplayName', 'DF Limit');
end
hold on
end