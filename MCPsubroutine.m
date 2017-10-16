function x = MCPsubroutine (u, mu, theta) %(x0,mu,theta) suffices
% solve min 0.5*norm(x-u)^2 + MCP(x) 
x = zeros(size(u,1),1);
s = sign(u); 


x2 = max(theta*mu, abs(u));
z0 = 0.5.*(x2-abs(u)).*(x2-abs(u))+(theta*mu^2/2);
z1 = 0.5.*(u.*u);
z2 = 0.5.*(theta*mu-abs(u)).*(theta*mu-abs(u))+(theta*mu^2/2);
w = min(theta*mu, max(0, theta*(abs(u))-mu)/(theta-1));
z3 = 0.5*(w-abs(u)).*(w-abs(u))+(mu.*abs(w)+(w.*w)/(2*theta));
z = min([z0, z1, z2, z3],[],2);
I = [z == z0];
x(I(:)) = x2(I(:));
% I = [z == z1];
% x(I(:)) = 0;
I = [z == z2];
x(I(:)) = mu*theta;
I = [z == z3];
x(I(:)) = w(I(:));
x = s.*x;
end
