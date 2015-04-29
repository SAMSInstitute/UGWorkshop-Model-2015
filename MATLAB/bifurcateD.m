% rect = [200 80 700 650]; %fix the window size and position
% set(0, 'defaultfigureposition',rect);
global N0
global D
 N0 = 4.5e-5;
% option=odeset('AbsTol',1e-11,'RelTol',1e-11);

tpts = linspace(7000,7500,501);


for D = 0.05:0.01:2.0
    Usol = ode45(@becks4Eq,[0,7500],[2.1e-7;4.56e-7;4.3e-5;4.5e-5]);
    u1 = deval(Usol, tpts, 1);
    u2 = deval(Usol, tpts, 2);
    %
    Figure1 = figure(1);
    set(Figure1, 'defaulttextinterpreter', 'latex')
%    set(gca, 'xlim', [0.05 2.0])
    hold on
    plot(D*ones(length(tpts),1), u1,'.','MarkerSize',2);
    title('Bifurcation Diagram for $R$ vs $D$');
    xlabel('Values for $D$');
    ylabel('Rod Species $R$');
%     subplot(2,1,2)
    Figure2 = figure(2);
    set(Figure2, 'defaulttextinterpreter', 'latex')
%    set(gca, 'xlim', [0.05 2.0])
    title('Bifurcation Diagram for $C$ vs $D$');
    xlabel('Values for $D$');
    ylabel('Cocci Species $C$');
    hold on
    plot(D*ones(length(tpts),1), u2,'.','MarkerSize',2);
    D
end
% ylabel('x','Rotation',0);
% set(gca, 'xlim', [0.5e-4 1.5e-4])
set(gca, 'xlim', [0.05 2.0])
hold off