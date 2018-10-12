function dudt = KSd_dirichlet(t0,u0,v)
% calculates the time derivative for spatially discretized modified KS eq:
% du/dt = -(u+c)du/dx-d^2u/dx^2-d^4u/dx^4
% with Dirichlet and Neumann boundary conditions
% u0 is the current state (row vectors); interior nodes only 
% c is the system free parameter

u0 = u0';

L = 16;
n = size(u0,2); dx = L/(n+1);
dudt = zeros(size(u0));

% -1/2*du^2/dx
dudt(:,2:end-1) = -(u0(:,3:end).^2-u0(:,1:end-2).^2)/4/dx;
dudt(:,1) = -u0(:,2).^2/4/dx;
dudt(:,end) = u0(:,end-1).^2/4/dx;

% -c*du/dx
% dudt(:,2:end-1) = dudt(:,2:end-1)-c*(u0(:,3:end)-u0(:,1:end-2))/2/dx;
% dudt(:,1) = dudt(:,1)-c*u0(:,2)/2/dx;
% dudt(:,end) = dudt(:,end)+c*u0(:,end-1)/2/dx;

% -d^2u/dx^2
dudt(:,2:end-1) = dudt(:,2:end-1)-(u0(:,1:end-2)-2*u0(:,2:end-1)...
                                    +u0(:,3:end))/dx^2;
dudt(:,1) = dudt(:,1)-(u0(:,2)-2*u0(:,1))/dx^2;
dudt(:,end) = dudt(:,end)-(u0(:,end-1)-2*u0(:,end))/dx^2;

% - d^4u/dx^4
dudt(:,3:end-2) = dudt(:,3:end-2)-v*(u0(:,1:end-4)-4*u0(:,2:end-3)+...
                    6*u0(:,3:end-2)-4*u0(:,4:end-1)+u0(:,5:end))/dx^4;
dudt(:,1) = dudt(:,1)-v*(7*u0(:,1)-4*u0(:,2)+u0(:,3))/dx^4;
dudt(:,end) = dudt(:,end)-v*(7*u0(:,end)-4*u0(:,end-1)+u0(:,end-2))/dx^4;
dudt(:,2) = dudt(:,2)-v*(-4*u0(:,1)+6*u0(:,2)-4*u0(:,3)+u0(:,4))/dx^4;
dudt(:,end-1) = dudt(:,end-1)-v*(-4*u0(:,end)+6*u0(:,end-1)-...
                                4*u0(:,end-2)+u0(:,end-3))/dx^4;
                            
end
