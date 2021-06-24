module SimpleCanonicalModel
    function m(alpha, gamma)
        return alpha./gamma
    end

    function p(alpha, beta, gamma, delta)
        return alpha.*beta./(gamma.*delta)
    end

    function m_var(alpha, gamma)
        return alpha./gamma
    end

    function p_var(alpha, beta, gamma, delta)
        return alpha.*beta.*(beta+gamma+delta)./(gamma.*delta.*(gamma+delta))
    end
end
