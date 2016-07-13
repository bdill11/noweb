function Aphilist2 = trans2d_Aphilist(T,Aphilist,order)
% function Aphilist2 = trans2d_Aphilist(T,Aphilist,order)
%
% Transforms Aphilist into Aphilist2 according to matrix T (2 x 2)
% Note: Order of colums is [val, d/dx1, d/dx2, d^2/dx1^2, d^2/dx1.dx2, d^2/dx2^2]
Aphilist2 = zeros(size(Aphilist));
Aphilist2(:,1) = Aphilist(:,1); % values unchanged
if order >= 1
    S = inv(T);
    Aphilist2(:,2:3) = Aphilist(:,2:3)*S; % chain rule for 1st derivatives
end % if
if order >= 2
    % chain rule for 2nd derivatives (affine transformation)
    Aphilist2(:,4) = Aphilist(:,4)*(S(1,1)^2)+ ...
           Aphilist(:,5)*(2*S(2,1)*S(1,1))+ ...
           Aphilist(:,6)*(S(2,1)^2);
    Aphilist2(:,5) = Aphilist(:,4)*(S(1,1)*S(1,2))+ ...
           Aphilist(:,5)*(S(1,1)*S(2,2)+S(1,2)*S(2,1))+ ...
           Aphilist(:,6)*(S(2,2)*S(2,1));
    Aphilist2(:,6) = Aphilist(:,4)*(S(1,2)^2)+ ...
           Aphilist(:,5)*(2*S(1,2)*S(2,2))+ ...
           Aphilist(:,6)*(S(2,2)^2);
end % if
end % function
