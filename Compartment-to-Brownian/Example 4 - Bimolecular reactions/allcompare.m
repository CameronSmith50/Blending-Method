close all


hold on

% 
% bar(phi_cube(1:100))
% %bar(phi_part(1:100),'r','facealpha',.0)
% plot((1.5:100.5),phi(1:100),'r')
% legend('compartment','analytical','location','northwest')


% figure
% hold on
plot(second_order)
plot(second_order1)
plot(second_order2)
plot(second_order_meso_cube_1COMP)

figure
plot((second_order-second_order_meso_cube_1COMP)/n0)