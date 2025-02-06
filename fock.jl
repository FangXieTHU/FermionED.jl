# Fock basis operations

module Fock

    using Printf

    const state_length = 64

    if state_length == 64
        bitmaplength = UInt64
    elseif state_length == 128
        bitmaplength = UInt128
    end

    const unsigned1::bitmaplength = 1
    const unsigned0::bitmaplength = 0

    mutable struct FermionFock
        sign::Int64
        state::bitmaplength
    end

    function init_state!()
        FermionFock(1, unsigned0)
    end

    function set_state!(fock_object, x)
        fock_object.sign = 1
        fock_object.state = x
        fock_object
    end

    function creation_operator!(fock_object, n)
        n_0based = n - 1
        state = fock_object.state
        # check if i-th site is occupied
        occ = (state >> n_0based) & unsigned1 == unsigned1
        if occ == false
            jw_coefficient = count_ones(state >> (n_0based + 1))
            if mod(jw_coefficient, 2) == 1
                fock_object.sign = fock_object.sign * (-1)
            end
            fock_object.state = state | (unsigned1 << n_0based)
        else
            fock_object.sign = 0
        end
        fock_object
    end

    function annihilation_operator!(fock_object, n)
        n_0based = n - 1
        state = fock_object.state
        # check if i-th site is occupied
        occ = (state >> n_0based) & unsigned1 == unsigned1
        if occ == true
            jw_coefficient = count_ones(state >> (n_0based + 1))
            if mod(jw_coefficient, 2) == 1
                fock_object.sign = fock_object.sign * (-1)
            end
            fock_object.state = xor(state, unsigned1 << n_0based)
        else
            fock_object.sign = 0
        end
        fock_object
    end

    # customize the print output of a state
    import Base.show
    function show(io::IO, fock_object::FermionFock)
        if fock_object.sign == 1
            s = bitstring(fock_object.state)
            output = @sprintf "+|%s⟩" s
        elseif fock_object.sign == -1
            s = bitstring(fock_object.state)
            output = @sprintf "-|%s⟩" s
        elseif fock_object.sign == 0
            s = bitstring(fock_object.state)
            output = "0"
        end
        print(output)
    end

end