% Function returns an extended vector, based on
% the type of boundary condition requested
% Vector is extended by m cells on each side


function U = apply_bc(Ui,bc,m)

switch bc
    case 'Periodic'
        U = [Ui(:,end-m+1:end) , Ui, Ui(:,1:m)];
    case 'Open'
        U = [repmat(Ui(:,1),1,m)   , Ui, repmat(Ui(:,end),1,m)];
end
