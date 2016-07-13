function [err_dphi,err_ddphi] = check_derivs(Aphifunc,x,d)
% function [err_dphi,err_ddphi] = check_derivs(Aphifunc,x,d)
%
% Returns errors in derivative test: err_dphi is the error vector for
% (phi(x+d)-phi(x-d)-2*grad phi(x)'*d)/norm(d),
% err_ddphi is the error vector for
% (grad phi(x+d)-grad phi(x-d)-2*Hess phi(x)*d)/norm(d).
%
% Assumes scalar element: order of rows:
% [phi(x), (d/dx1)phi(x), (d/dx2)phi(x), (d^2/dx1^2)phi(x), ...
% (d^2/dx1.dx2)phi(x), (d^2/dx2^2)phi(x)]
Aphivalx   = Aphifunc(x,2);
Aphivalxpd = Aphifunc(x+d,2);
Aphivalxmd = Aphifunc(x-d,2);
phixpd = Aphivalxpd(:,1);
phixmd = Aphivalxmd(:,1);
dphix  = Aphivalx(:,2:3);
err_dphi = (phixpd-phixmd-2*dphix*d)/norm(d);
dphixpd = Aphivalxpd(:,2:3);
dphixmd = Aphivalxmd(:,2:3);
ddphix  = Aphivalx(:,4:6);
err_ddphi = (dphixpd-dphixmd-2*(d(1)*ddphix(:,1:2)+d(2)*ddphix(:,2:3)))/norm(d);
