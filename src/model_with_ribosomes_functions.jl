module ModelWithRibosomes
    export mrna, p, r, mrna_var, r_var, p_var, m_r_covar, p_r_covar, m_p_covar
    function mrna(alpha, gamma, T, t)
        return alpha/gamma*(1-exp(-gamma * t)/(2-exp(-gamma*T)))
    end

    function p(alpha, beta, gamma, rho, T, t)
        return alpha*beta*(2+4*T*gamma - (t^2 + 2*t*T + 3*T^2)*gamma^2 + 2*exp((T-t)*gamma)*(1+(t+T)*gamma)+2*exp(T*gamma)*(-2+gamma*(t^2*gamma+T*(-2+2*t*gamma+3*T*gamma))))*rho / (2 * (-1 + 2*exp(T*gamma))*gamma^3)
    end

    function r(rho, T, t)
        return rho*(t+T)
    end

    function mrna_var(alpha, gamma, T, t)
        return alpha*(exp(gamma*(T-t)) - 2*exp(gamma*T) + 1)/(gamma*(1-2*exp(gamma*T)))
    end

    function r_var(rho, T, t)
        return rho*(t+T)
    end

    function p_var(alpha, beta, gamma, rho, T, t)
        return (alpha*beta*rho*exp(-gamma*t))/(18*gamma^5*(1-2*exp(gamma*T))^2)*(exp(gamma*t)*(alpha*beta*(2*gamma*(gamma^2*(9*t^2*T+3*t^3+9*t*T^2+7*T^3)-4*gamma*T*(3*t+4*T)-6*t+8*T)+25)-2*beta*rho*(gamma*(3*gamma*t^2*(3-2*gamma*t)+6*gamma*T^2*(7-3*gamma*t)-18*T*(gamma*t-2)*(gamma*t+1)+4*gamma^2*T^3)+18)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-4))-2))-4*exp(gamma*(t+T))*(beta*rho*(gamma*(6*gamma*t^2*(2*gamma*t-3)+9*gamma*T^2*(4*gamma*t-7)+36*T*(gamma*t-2)*(gamma*t+1)+10*gamma^2*T^3)-54)+2*alpha*beta*(3*gamma^3*t^3+9*gamma^3*t^2*T+3*gamma*t*(gamma*T*(3*gamma*T-4)-3)+gamma*T*(gamma*T*(7*gamma*T-13)+3)+13)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-3))-2))+4*exp(gamma*(t+2*T))*(alpha*beta*(2*gamma*(gamma^2*(9*t^2*T+3*t^3+9*t*T^2+7*T^3)-2*gamma*T*(6*t+5*T)-12*t+T)+31)+2*beta*rho*(gamma*(3*gamma*t^2*(2*gamma*t-3)+3*gamma*T^2*(6*gamma*t-7)+18*T*(gamma*t-2)*(gamma*t+1)+14*gamma^2*T^3)-36)+9*gamma^2*(gamma*(gamma*t^2+T*(2*gamma*t+3*gamma*T-2))-2))+12*exp(2*gamma*T)*(18*beta*rho+3*gamma*(gamma+beta*gamma*rho*(3*t^2+6*t*T+2*T^2)+gamma^2*(t+T)*(beta*rho*t*(t+2*T)+1)+6*beta*rho*(t+T))+alpha*beta*(3*gamma^2*(t+T)^2-4*gamma*T-10))-6*exp(gamma*T) *(18*beta*rho+3*gamma *(gamma+beta*gamma*rho*(3*t^2+6*t*T+2*T^2)+gamma^2*(t+T)*(beta*rho* t*(t+2*T)+1)+6*beta*rho*(t+T))+alpha*beta*(3*gamma^2*(t+T)^2-4*gamma*T-8))+9*alpha*beta*exp(-(gamma*(t-2*T)))*(2*gamma*(t+T)+3))
    end

    function m_r_covar()
        return 0
    end

    function p_r_covar(alpha, beta, gamma, rho, T, t)
       return alpha*beta*rho/(6*gamma^3*(2*exp(gamma*T)-1))*(-3*gamma^2*(t+T)^2+2*exp(gamma*T)*(3*gamma^2*(t+T)^2 - 4*gamma*T - 4) + 6*exp(gamma*(T-t))*(gamma*(t+T)+1)+4*gamma*T+2) 
    end

    function m_p_covar(alpha, beta, gamma, rho, T, t)
        return exp(-t*gamma)*alpha*beta*rho*(-exp(gamma*T)*t*(t+2*T)*gamma^2+2*exp(gamma*t)*(gamma*(t+T)-1)*(2*exp(gamma*T)-1))/(2*(2*exp(T*gamma)-1)*gamma^3)
    end

    function p_time_av(alpha, beta, gamma, rho, T)
        return alpha*beta*(-12 + T*gamma *(-6+T*gamma*(12-13*T*gamma))+2*exp(gamma*T)*(6+gamma*T*(-3+gamma*T*(-6 + 13*gamma*T))))*rho/(6*(-1+2*exp(gamma*T))*T*gamma^4)
    end

    function p_var_time_av(alpha, beta, gamma, rho, T)
       (alpha*beta*rho*(34*alpha*beta*gamma^6*rho*T^6*(1-2*exp(gamma*T))^2+5*beta*gamma^6*T^5*(2*exp(gamma*T)-1)*(-61*alpha-50*rho+122*(alpha+2*rho)*exp(gamma*T))-10*gamma^4*T^4*(2*exp(gamma*T)-1)*(-44*alpha*beta*gamma+60*alpha*beta*rho+39*gamma^3+2*exp(gamma*T)*(4*alpha*beta*(8*gamma+3*rho)+66*beta*gamma*rho-39*gamma^3))-10*gamma^3*T^3*(alpha*beta*(-82*gamma-168*rho+(4*gamma+15*rho)*exp(2*gamma*T)+6*(25*gamma+66*rho)*exp(gamma*T))+36*gamma*(2*exp(gamma*T)-1)*(beta*rho*(2*exp(gamma*T)+9)+gamma^2*(exp(gamma*T)-1)))+10*gamma^2*T^2*(alpha*beta*(55*gamma+162*rho-4*(53*gamma+162*rho)*exp(gamma*T)+(157*gamma+297*rho)*exp(2*gamma*T))-18*gamma*(2*exp(gamma*T)-1)*(34*beta*rho+gamma^2+exp(gamma*T)*(gamma^2-10*beta*rho)))+15*gamma*T*(exp(gamma*T)-1)*(alpha*beta*(20*gamma+15*rho+(159*rho-20*gamma)*exp(gamma*T))+24*gamma*(2*exp(gamma*T)-1)*(12*beta*rho+gamma^2))-720*alpha*beta*rho*(exp(gamma*T)-1)^2))/(180*gamma^8*(T-2*T*exp(gamma*T))^2)
    end
end
