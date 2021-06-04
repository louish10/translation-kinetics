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
       -(alpha*beta*rho*exp(gamma*T)*(-8*beta*rho*(alpha*(gamma^3*T^3*(gamma*T*(2*gamma*T*(431*gamma*T-585)+675)+1665)+270*gamma*T+180)-5*gamma^2*T*(gamma*T*(gamma*T*(gamma*T*(43*gamma*T-33)+144)+351)+324))+gamma*T*sinh(gamma*T)*(5*gamma^2*(-72*gamma+alpha*beta*T*(gamma*T*(3*gamma*T*(56-61*gamma*T)+172)-204)+18*gamma^2*T*(gamma*T*(4-13*gamma*T)+6))+6*beta*rho*(alpha*(gamma*T*(gamma*T*(2*gamma*T*(gamma*T*(431*gamma*T-260)-480)+695)+225)-45)+5*gamma*(gamma*T*(gamma*T*(gamma*T*(88-73*gamma*T)+156)+84)-144)))+20*gamma^2*T*(alpha*beta*(gamma*T*(gamma*T*(gamma*T*(61*gamma*T-76)+75)+106)-30)+3*gamma^2*(gamma*T*(2*gamma*T*(13*gamma*T-9)+3)+18))+5*cosh(gamma*T)*(alpha*beta*gamma^2*T*(120-gamma*T*(gamma*T*(gamma*T*(305*gamma*T-344)+156)+424))+2*alpha*beta*rho*(gamma*T*(gamma*T*(gamma*T*(2*gamma*T*(gamma*T*(431*gamma*T-468)+72)+1089)-45)+216)+144)-2*beta*gamma^2*rho*T*(gamma*T*(gamma*T*(gamma*T*(269*gamma*T-264)+180)+972)+1296)-6*gamma^4*T*(gamma*T*(gamma*T*(65*gamma*T-36)-6)+36))))/(180*gamma^8*(T-2*T*exp(gamma*T))^2) 
    end
end
