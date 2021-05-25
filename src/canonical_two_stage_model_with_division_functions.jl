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
end
