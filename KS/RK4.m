function u = RK4(dudt,u0,t0,dt)

    k1 = dudt(t0,u0);
    k2 = dudt(t0+dt/2,u0+dt*k1/2);
    k3 = dudt(t0+dt/2,u0+dt*k2/2);
    k4 = dudt(t0+dt,u0+dt*k3);
    
    u = u0 + 1/6*(k1+2*k2+2*k3+k4)*dt;

end
