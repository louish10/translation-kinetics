module CanonicalTwoStageModel
    export mrna, approximate_protein_variance, approximate_protein_mean, approximate_covariance

    function mrna(alpha, gamma, T, t)
        return alpha/gamma*(1-exp(-gamma * t)/(2-exp(-gamma*T)))
    end

    function approximate_covariance(alpha, beta, gamma, T, t)
        D = 1/(4-exp(-gamma*T))*(3/gamma + T*exp(-gamma * T)/(2-exp(-gamma*T)))
        return beta*alpha/gamma*(1/gamma - t*exp(-gamma*t)/(2-exp(-gamma*T)) - exp(-gamma * t)*D)
    end

    function approximate_protein_mean(alpha, beta, gamma, T, t)
        A = alpha/gamma*(1/(2-exp(-gamma * T)))
        c = beta*alpha*T/gamma - 2*beta*A/gamma + beta*A/gamma*exp(-gamma * T)
        return beta*alpha/gamma * (t+T+1/gamma*(exp(-gamma*t) + exp(-gamma * T) - 2)/(2-exp(-gamma*T)))
    end

    function approximate_protein_variance(alpha, beta, gamma, T, t)
        c_double_prime = (approximate_F(alpha, beta, gamma, T, T)- 4 * approximate_F(alpha, beta, gamma, T, 0) + approximate_protein_mean(alpha, beta, gamma, T, T))/3
        return approximate_F(alpha, beta, gamma, T, t) + c_double_prime
    end

    function approximate_F(alpha, beta, gamma, T, t)
        A = alpha/gamma - mrna(alpha, gamma, T, 0)
        D = 1/(4-exp(-gamma*T))*(3/gamma + T*exp(-gamma * T)/(2-exp(-gamma*T)))
        return beta*alpha/gamma * ((2*beta/gamma+1)*t + exp(-gamma*t)/(2-exp(-gamma*T))*(2*beta*(gamma*t + 1)+gamma)/(gamma^2) + 2*beta*exp(-gamma*t)/gamma*D)
    end

    function m_var_time_av(alpha, gamma, T)
        (alpha*(alpha*(exp(gamma*T)-1)*(gamma*T+exp(gamma*T)*(gamma*T-2)+2)+2*gamma^2*T*(2*exp(gamma*T)-1)*(gamma*(-T)+exp(gamma*T)*(2*gamma*T-1)+1)))/(2*gamma^4*(T-2*T*exp(gamma*T))^2)
    end

    function p_time_av(alpha, beta, gamma, T)
        return alpha*beta*(1+1/(1-2*exp(gamma*T))+gamma*T*(-2+3*T*gamma))/(2*T*gamma^3)
    end

    function p_var_time_av(alpha, beta, gamma, T)
        return (alpha  *beta  *(alpha  *beta*(gamma *T *(gamma^3*T^3+12 *gamma*T+18)-4*exp(gamma *T) *(gamma *T*(gamma^3*T^3+3 *gamma*T+18)-6)+2*exp(2 *gamma *T)*(gamma *T *(2 *gamma^3*T^3-12* gamma *T+27)-6)-12)+2*gamma ^2*T *(2*exp(gamma* T)-1)*(-2 *beta *(gamma *T*(3 *gamma*T+4)+12)+4 *beta *exp(gamma *T)*(gamma *T*(5*gamma *T-4)+6)+3 *gamma  *(gamma *T*(2-3 *gamma *T)+exp(gamma *T)*(2 *gamma *T*(3*gamma *T-2)+2)-2))))/(12*gamma^6 *(T-2*T*exp(gamma *T))^2)
    end
end
