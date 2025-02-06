# recursively generate the hilbert space with provided good quantum numebers

module HilbertSpace
    
    import Base.+
    import Base.==

    const state_length = 64

    if state_length == 64
        bitmaplength = UInt64
    elseif state_length == 128
        bitmaplength = UInt128
    end
    
    unsigned1::bitmaplength = 1

    mutable struct GoodQuantumNumbers
        L1::Int64
        L2::Int64
        k1::Int64
        k2::Int64
        singlebody_quantum_numbers::Matrix{Int64}
    end

    function init_good_quantum_numbers!(L1, L2, good_quantum_numbers)
        gqn = GoodQuantumNumbers(L1, L2, 0, 0, good_quantum_numbers)
        return gqn
    end

    function +(Q1::GoodQuantumNumbers, dQ::Vector{Int64})
        k1_tot = mod(Q1.k1 + dQ[1], Q1.L1)
        k2_tot = mod(Q1.k2 + dQ[2], Q1.L2)
        Q2 = GoodQuantumNumbers(Q1.L1, Q1.L2, k1_tot, k2_tot, Q1.singlebody_quantum_numbers)
        return Q2
    end

    function ==(Q1::GoodQuantumNumbers, target::Vector{Int64})
        if mod(Q1.k1 - target[1], Q1.L1) == 0 && mod(Q1.k2 - target[2], Q1.L2) == 0
            return true
        else
            return false
        end
    end

    function reset_good_quantum_number!(Q1::GoodQuantumNumbers)
        Q1.k1 = 0
        Q1.k2 = 0
    end

    function add_one_particle!(N_size, sec)
        sec_len = length(sec)
        for i in 1:sec_len
            sec[i] = sec[i] + (unsigned1 << (N_size - 1))
        end
    end
    
    function hilbert_space(N_particle, gqn, goal)
        reset_good_quantum_number!(gqn)
        N_size = size(gqn.singlebody_quantum_numbers)[1]
        hilbert_space_recursive(N_size, N_particle, gqn, goal)
    end

    function hilbert_space_recursive(N_size, N_particle, qgn, goal)
        if N_size < N_particle || N_particle < 0
            return Array{bitmaplength}(undef, 0)
        elseif N_particle == 0
            is_match = qgn == goal
            if is_match
                return Array{bitmaplength}([0])
            else
                return Array{bitmaplength}(undef, 0)
            end
        else
            sec1 = hilbert_space_recursive(N_size-1, N_particle, qgn, goal)
            sec2 = hilbert_space_recursive(N_size-1, N_particle-1, qgn+qgn.singlebody_quantum_numbers[N_size,:], goal)
            add_one_particle!(N_size, sec2)
            append!(sec1, sec2)
        end
    end

    # generate the dictionary for search the fock state index
    function generate_hashtable(hs_array)
        l = length(hs_array)
        hs_dict = Dict(hs_array[i]=>i for i in 1:l)
        hs_dict
    end

    function find_target_state(hs_dict, fock_state_bin)
        hs_index = get(hs_dict, fock_state_bin, -1)
        hs_index
    end
end
