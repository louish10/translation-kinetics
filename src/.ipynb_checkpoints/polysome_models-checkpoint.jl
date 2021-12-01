module PolysomeModelA
   export m1, m2, r, r1

    function m1(k1, k2)
        return k1/k2
    end

    function m2(k1p, k2p)
        return k1p/k2p
    end

    function r(k7, k8)
        return k7/k8
    end
    
    function r1(k1, k2, k3, k4, k5, k6, k7, k8)
        return k1*k3*k7/(k2*k4*k8 + k2*k5*k8)
    end

    function r2(k1p, k2p, k3p, k4p, k5p, k6p, k7, k8)
        return k1p*k3p*k7/(k2p*k4p*k8 + k2p*k5p*k8)
    end

    function p1(k1, k2, k3, k4, k5, k6, k7, k8)
        return k1*k3*k5*k7/(k2*k4*k6*k8+k2*k5*k6*k8)
    end

    function p2(k1p, k2p, k3p, k4p, k5p, k6p, k7, k8)
        return k1p*k3p*k5p*k7/(k2p*k4p*k6p*k8+k2p*k5p*k6p*k8)
    end

    function m1_var(k1, k2)
        return k1/k2
    end
    
    function r_var(k7,k8)
        return k7/k8
    end

    function r1_var(k1,k2,k3,k4,k5,k7,k8)
        return k1*k3*k7/(k2*k4*k8+k2*k5*k8)
    end

    function p1_var(k1,k2,k3,k4,k5,k6,k7,k8)
    return (k1*k3*k5*k7*(k1*k3*k6*(k5+k6)*k8+k2^2*(k4+k5+k6)*k8*(k6+k8)+k2*k6*(k4+k5+k6)*k8*(k6+k8)+k2*k3*(k5+k6)*(k6*k7+(k1+k7)*k8)))/(k2*(k4+ k5)*k6*k8*(k1*k3*k6^2*k8+k2^2*(k4+k5+k6)*k8*(k6+k8)+k2*k6*(k3*k6*k7+k3*(k1+k7)*k8+(k4+k5+k6)*k8*(k6+k8))))
    end
end

module PolysomeModelB
    function m1(k1,k2)
        k1/k2
    end

    function r(k7,k8)
        k7/k8
    end

    function r1(k1, k2, k3, k4, k5, k6, k7, k8)
        return k1*k3*k7/(k2*k4*k8 + k2*k5*k8)
    end
    
    function p1(k1, k2, k3, k4, k5, k6, k7, k8)
        return k1*k3*k5*k7/(k2*k4*k6*k8+k2*k5*k6*k8)
    end

    function m1_var(k1,k2)
        k1/k2
    end

    function r_var(k7,k8)
        k7/k8
    end

    function r1_var(k1,k2,k3,k4,k5,k7,k8)
        return k1*k3*k7/(k2*k4*k8+k2*k5*k8)
    end

    function p1_var(k1,k2,k3,k4,k5,k6,k7,k8)
        return (k1*k3*k5*k7*(k1*k3*k6*(k5+k6)*k8+k2^2*(k4+k5+k6)*k8*(k6+k8)+k2*k6*(k4+k5+k6)*k8*(k6+k8)+k2*k3*(k5+k6)*(k6*k7+(k1+k7)*k8)))/(k2*(k4+ k5)*k6*k8*(k1*k3*k6^2*k8+k2^2*(k4+k5+k6)*k8*(k6+k8)+k2*k6*(k3*k6*k7+k3*(k1+k7)*k8+(k4+k5+k6)*k8*(k6+k8))))
    end

end

module PolysomeModelC
    export m1, r, p1, m1_var, r_var, p1_var

    function m1(k1, k2)
        return k1/k2
    end

    
    function r(k7, k8)
        return k7/k8
    end

    function p1(k1, k2, k3, k6, k7, k8)
        return k1*k3*k7/(k2*k6*k8)
    end

    function m1_var(k1, k2)
        return k1/k2
    end
    
    function r_var(k7,k8)
        return k7/k8
    end

    function p1_var(k1,k2,k3,k6,k7,k8)
        return k1*k3*k7*(k2*k3*k6*k7+(k2+k6)*(k1*k3+k2*k6)*k8+k2*k3*k7*k8+k2*(k2+k6)*k8^2)/(k2^2*k6*(k2+k6)*k8^2*(k6+k8))
    end
end