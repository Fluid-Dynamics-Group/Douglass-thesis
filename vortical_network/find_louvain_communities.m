function [Ci] = find_louvain_communities(gamma,A)
    % find_louvain_communities
    % Function to find communities in vortical (or any) networks

    % Using Modularity max.
    % [Ci, Q] = modularity_dir(A,gamma);

    % Using Louvain algorithm
    Ci = 1:length(A);               % initial community affiliations                   
    Q0 = -1;                        % initialize modularity values
    Q = 0;                 
    while Q-Q0>1e-5                 % while modularity increases
        Q0 = Q;                     % perform community detection
        [Ci, Q] = community_louvain(A,gamma,Ci);
    end
end
