%     Copyright 2014 Hugo Bowne-Anderson
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in this script, we write a function that solves the Gierer_Meinhardt 
%system for parameters D, mu and sigma. We use MATLAB's PDEPE to do so.
%The function also plots the solution 1) as a surface & 2) as an animation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = GM_solve(P,tmax,n,L,fig_numbers)
    %we solve the PDE for parameters
    %P(1) = D in PDE;
    %P(2) = mu in PDE;
    %P(3) = sigma in PDE;
    %tmax = amount of time that we run the system;
    %n = number of time points that we integrate over
    %L = length (spatial) of system
    
    GMfunctions = GMfuns; %here we call the necessary functions to solve the PDE:
                     %the initial conditions function, the boundary conditions
                     %function and the PDE function itself.
    rng(404); %set the seed for random number generation
    t = linspace(0,tmax,n); %time vector
    x = linspace(0,L,200); %the mesh on which we solve
    %%PDEPE solver
    sol = pdepe(0,@GMfunctions.fun3,@GMfunctions.fun1,@GMfunctions.fun2,x,t,[],P);
    u1 = sol(:,:,1); %activator solution
    u2 = sol(:,:,2); %inhibitor solution



    %% Plotting in 3 DIMENSIONS

    figure(fig_numbers(1))
    surf(x,t,u1,'edgecolor','none');
    xlabel('Position','fontsize',20,'fontweight','b','fontname','arial')
    ylabel('Time','fontsize',20,'fontweight','b','fontname','arial')
    zlabel('[Activator]','fontsize',20,'fontweight','b','fontname','arial')
    axis([0 L 0 tmax 0 max(max(u1))])
    set(gcf(), 'Renderer', 'painters')
    set(gca,'FontSize',18,'fontweight','b','fontname','arial')

    figure(fig_numbers(2))
    surf(x,t,u2,'edgecolor','none');
    xlabel('Position','fontsize',20,'fontweight','b','fontname','arial')
    ylabel('Time','fontsize',20,'fontweight','b','fontname','arial')
    zlabel('[Inhibitor]','fontsize',20,'fontweight','b','fontname','arial')
    axis([0 L 0 tmax 0 max(max(u2))])
    set(gcf(), 'Renderer', 'painters')
    set(gca,'FontSize',18,'fontweight','b','fontname','arial')

    %%

    figure(fig_numbers(3))
    for n = 1:length(t)
        set(gca, 'FontSize', 18, 'LineWidth', 1); %<- Set properties
        plot( x , sol(n,:,1), 'LineWidth',3);
        hold on
        plot( x , sol(n,:,2), 'r', 'LineWidth',3);
        hold off
        legend('activator', 'inhibitor', 'Location', 'SouthEast');
        title(strcat('Gierer-Meinhardt patterns t =' , sprintf(' %d ', ceil(t(n)))));
        axis([0 L 0 max(max(max(sol(:,:,:))))+0.1])
        M(n) = getframe;
    end

    % play as smooth movie 1 time at 5 frames per second
    numtimes=1;
    fps=5;
    movie(M,numtimes,fps)
end