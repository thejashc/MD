function [ ] = visualize(r)

global box

% visualize the positions of the particles
plot(r(:,1), r(:,2), 'o','Markerfacecolor','r','MarkerEdgecolor','r','MarkerSize',5);
%scatter(r(:,1), r(:,2))
xlim([-box/2 box/2])
ylim([-box/2 box/2])

grid on;

end

