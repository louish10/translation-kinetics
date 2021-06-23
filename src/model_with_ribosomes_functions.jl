module ModelWithRibosomes
    export mrna, p, r, mrna_var, r_var, p_var, m_r_covar, p_r_covar, m_p_covar, limit_of_differences
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

    function m_time_av(alpha, gamma, T)
        alpha*(gamma*T*(1-2*exp(gamma*T)) - 1 + exp(gamma*T))/(T*gamma^2*(1-2*exp(gamma*T)))
    end

    function m_var_time_av(alpha, gamma, T)
        return (alpha  *(alpha  *(exp(gamma *T)-1)*(gamma *T+exp(gamma *T)*(gamma *T-2)+2)+2*gamma ^2*T *(2*exp(gamma *T)-1)*(gamma *(-T)+exp(gamma *T)*(2 *gamma *T-1)+1)))/(2*gamma ^4 *(T-2*T*exp(gamma *T))^2)
    end

    function p_var_time_av(alpha, beta, gamma, rho, T)
        return (alpha  *beta *rho  *(beta *(5 *rho  *(alpha  *(4*exp(gamma *T)-1)*(gamma *T*(gamma *T*(13 *gamma *T-12)+6)-2*exp(gamma *T)*(gamma *T*(gamma *T*(13*gamma *T-6)-3)+6)+12)^2+2 *gamma ^2*T *(gamma *T*(gamma *T*(gamma *T*(61 *gamma *T-138)+408)+576)-16*exp(3 *gamma *T)*(gamma *T*(gamma *T*(gamma *T*(61 *gamma *T-90)-6)+135)+162)-2*exp(gamma *T)*(gamma *T*(gamma *T*(gamma *T*(208 *gamma *T-291)+1218)+1944)+972)+4*exp(2 *gamma *T)*(gamma *T*(gamma *T*(gamma *T*(269 *gamma *T-333)+798)+1638)+1080)+216))-alpha  *gamma *T*exp((3 *gamma *T)/(2))*(5 *sinh((gamma *T)/2)+3 *cosh((gamma* T)/2)) *((5 *gamma ^2*T*(gamma *T*(3 *gamma *T*(61 *gamma *T-56)-172)+204)+3 *rho*(gamma *T*(gamma *T*(gamma *T*(gamma *T*(879 *gamma *T-520)-1140)+390)+450)+390)) *sinh(gamma *T)-12 *rho *(gamma *T*(gamma *T*(gamma *T*(gamma *T*(293 *gamma *T-390)+265)+720)+270)+180)+5*(432 *rho +gamma *(3 *rho *T*(gamma *T*(gamma *T*(gamma *T*(293 *gamma *T-312)+52)+414)+138)+gamma *T*(gamma *T*(gamma *T*(305 *gamma *T-344)+156)+424)-120))*cosh(gamma* T)-20 *gamma *(gamma *T*(gamma *T*(gamma *T*(61 *gamma *T-76)+75)+106)-30)))-30 *gamma ^4*T *(-6*exp(gamma *T)+8*exp(2 *gamma *T)+1) *(gamma *T*(gamma *T*(12-13 *gamma *T)-6)+2*exp(gamma *T)*(gamma *T*(gamma *T*(13 *gamma *T-6)-3)+6)-12)))/(180 *gamma ^8 *(1-4*exp(gamma *T)) *(T-2*T*exp(gamma *T))^2)
    end

    function limit_of_differences(alpha, beta, gamma, T)
        return (beta *(alpha  *(4*exp(gamma *T)-1) *(gamma *T *(gamma *T*(gamma *T *(gamma ^3*T^3+780*gamma *T+4290)+8100)-900)-4*exp(gamma *T) *(gamma *T *(gamma *T*(gamma *T *(gamma ^3*T^3+555*gamma *T+1530)+7290)+2160)-1440)+2*exp(2 *gamma *T) *(gamma *T*(gamma *T *(gamma *T *(2*gamma ^3*T^3+660 *gamma *T-3945)+6750)+4770)-1440)-2880)+20 *gamma ^2*T *(gamma *T*(gamma *T*(gamma *T*(13 *gamma *T+60)-6)-1152)+exp(gamma *T)*(2 *gamma *T*(gamma *T*(31 *gamma *T*(3-4 *gamma *T)+87)+3888)+3888)-8*exp(3 *gamma *T)*(gamma *T*(gamma *T*(2 *gamma *T*(13 *gamma*T-90)+591)-540)-648)+4*exp(2 *gamma *T)*(gamma *T*(gamma *T*(gamma *T*(137 *gamma *T-333)+510)-3276)-2160)-432))-90*gamma ^4*T^2 *(-6*exp(gamma *T)+8*exp(2*gamma *T)+1) *(gamma *(-T)*(gamma *T*(gamma *T+6)+6)+2*exp(gamma *T)*(gamma *T*(gamma *T*(gamma *T-6)+15)-12)+24))/(135 *gamma ^2*T^2*(beta *(alpha  *(gamma *(-T) *(gamma ^3*T^3+12 *gamma*T+18)+8*exp(3 *gamma *T) *(gamma *T*(2 *gamma ^3*T^3-12 *gamma*T+27)-6)+4*exp(gamma *T)*(gamma *T *(2 *gamma ^3*T^3+15*gamma *T+36)-18)-2*exp(2 *gamma*T) *(gamma *T *(10 *gamma ^3*T^3+12*gamma *T+171)-54)+12)+4*gamma ^2*T *(2*exp(gamma *T)-1)*(gamma *T*(5 *gamma *T-8)+exp(2 *gamma *T)*(40 *gamma *T*(gamma *T-2)+84)-2*exp(gamma* T)*(gamma *T*(15 *gamma *T-17)+57)+30))+6 *gamma ^3*T*(-6*exp(gamma *T)+8*exp(2 *gamma *T)+1) *(gamma *T*(2-3 *gamma *T)+exp(gamma *T)*(2 *gamma *T*(3 *gamma*T-2)+2)-2)))
    end
end
