% rect = [200 80 700 650]; %fix the window size and position
% set(0, 'defaultfigureposition',rect);
global N0

% option=odeset('AbsTol',1e-11,'RelTol',1e-11);

tpts = linspace(7000,7500,501);

% for N0 = 1.5e-7:1.0e-8:5.0e-7  %Sherli's data
for N0 = 1.5e-5:1.0e-7:4.5e-5
% for N0 = 0.5e-4:1.0e-6:1.5e-4
% for N0=0.5e-4:1.0e-4:2.5e-4
%     [t,U] = ode45('becks4Eq',[0,7500],[2.1e-7;4.56e-7;4.3e-5;N0]);
    Usol = ode45(@becks4Eq,[0,7500],[2.1e-7;4.56e-7;4.3e-5;N0]);
%     Usol = ode45(@becks4Eq,[0,7500],[8.6703e-8;6.8236e-7;5.5372e-8;N0]);
%    Sherli's data above
    u1 = deval(Usol, tpts, 1);
    u2 = deval(Usol, tpts, 2);
    %
%     subplot(2,1,1)
    Figure1 = figure(1);
    set(Figure1, 'defaulttextinterpreter', 'latex')
    hold on
    plot(N0*ones(length(tpts),1), u1,'.','MarkerSize',2);
    title('Bifurcation Diagram for $R$ vs $N_0$');
    xlabel('Values for $N_0$');
    ylabel('Rod Species $R$');
%     subplot(2,1,2)
    Figure2 = figure(2);
    set(Figure2, 'defaulttextinterpreter', 'latex')
    title('Bifurcation Diagram for $C$ vs $N_0$');
    xlabel('Values for $N_0$');
    ylabel('Cocci Species $C$');
    hold on
    plot(N0*ones(length(tpts),1), u2,'.','MarkerSize',2);
    % axis(ax);
%     xlim([1.5e-5 4.5e-5]);
    N0
end
set(gca, 'xlim', [1.5e-5 4.5e-5])
hold off