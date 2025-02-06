# construct the interacting many body Hamiltonian

module ManyBodyHamiltonian

    using SparseArrays

    include("fock.jl")
    include("hilbert_space.jl")

    # determine the length of csr matrices elements
    function hamiltonian_data_length(hilbert_space_array, ts, t_c1dag, t_c2, Vs, V_c1dag, V_c2dag, V_c3, V_c4, kinetic)
        data_length = 0
        num_t_terms = length(ts)
        num_V_terms = length(Vs)
        dim = length(hilbert_space_array)

        # declare a FermionFock object and set it to defaul value
        fock_state = Fock.init_state!()

        for i in 1:dim
            x = hilbert_space_array[i]

            for V_term in 1:num_V_terms
                Fock.set_state!(fock_state, x)
                pos4 = V_c4[V_term]
                pos3 = V_c3[V_term]
                pos2 = V_c2dag[V_term]
                pos1 = V_c1dag[V_term]

                Fock.annihilation_operator!(fock_state, pos4)
                if fock_state.sign == 0
                    continue
                end
                Fock.annihilation_operator!(fock_state, pos3)
                if fock_state.sign == 0
                    continue
                end
                Fock.creation_operator!(fock_state, pos2)
                if fock_state.sign == 0
                    continue
                end
                Fock.creation_operator!(fock_state, pos1)
                if fock_state.sign == 0
                    continue
                end
                data_length += 1
            end
            
            if kinetic
                for Ek_term in 1:num_t_terms
                    Fock.set_state!(fock_state, x)
                    pos2 = t_c2[Ek_term]
                    pos1 = t_c1dag[Ek_term]

                    Fock.annihilation_operator!(fock_state, pos2)
                    if fock_state.sign == 0
                        continue
                    end
                    Fock.creation_operator!(fock_state, pos1)
                    if fock_state.sign == 0
                        continue
                    end
                    data_length += 1
                end
            end
        end
        data_length
    end

    function hamiltonian_sparse_data(data_length, hilbert_space_array, hilbert_space_dict, ts, t_c1dag, t_c2, Vs, V_c1dag, V_c2dag, V_c3, V_c4, kinetic)

        num_t_terms = length(ts)
        num_V_terms = length(Vs)
        dim = length(hilbert_space_array)

        # indptr = zeros(Int64, dim+1)
        # indptr[1] = 1
        col_indices = zeros(Int64, data_length)
        row_indices = zeros(Int64, data_length)
        data = zeros(ComplexF64, data_length)

        data_index = 1

        # declare a FermionFock object and set it to defaul value
        fock_state = Fock.init_state!()

        for i in 1:dim
            x = hilbert_space_array[i]

            for V_term in 1:num_V_terms
                Fock.set_state!(fock_state, x)
                pos4 = V_c4[V_term]
                pos3 = V_c3[V_term]
                pos2 = V_c2dag[V_term]
                pos1 = V_c1dag[V_term]
                V = Vs[V_term]
                
                Fock.annihilation_operator!(fock_state, pos4)
                if fock_state.sign == 0
                    continue
                end
                Fock.annihilation_operator!(fock_state, pos3)
                if fock_state.sign == 0
                    continue
                end
                Fock.creation_operator!(fock_state, pos2)
                if fock_state.sign == 0
                    continue
                end
                Fock.creation_operator!(fock_state, pos1)
                if fock_state.sign == 0
                    continue
                end

                j = HilbertSpace.find_target_state(hilbert_space_dict, fock_state.state)
                if j != -1
                    row_indices[data_index] = j
                    col_indices[data_index] = i
                    data[data_index] = V * fock_state.sign
                    data_index += 1
                end

            end
            
            if kinetic
                for Ek_term in 1:num_t_terms
                    Fock.set_state!(fock_state, x)
                    pos2 = t_c2[Ek_term]
                    pos1 = t_c1dag[Ek_term]
                    t = ts[Ek_term]

                    Fock.annihilation_operator!(fock_state, pos2)
                    if fock_state.sign == 0
                        continue
                    end
                    Fock.creation_operator!(fock_state, pos1)
                    if fock_state.sign == 0
                        continue
                    end

                    j = HilbertSpace.find_target_state(hilbert_space_dict, fock_state.state)
                    if j != -1
                        row_indices[data_index] = j
                        col_indices[data_index] = i
                        data[data_index] = t * fock_state.sign
                        data_index += 1
                    end
                end
            end
            # indptr[i + 1] = data_index
        end
        return row_indices, col_indices, data
    end

    function hamiltonian(hilbert_space_array, hilbert_space_dict, ts, t_c1dag, t_c2, Vs, V_c1dag, V_c2dag, V_c3, V_c4, kinetic)
        data_length = hamiltonian_data_length(hilbert_space_array, ts, t_c1dag, t_c2, Vs, V_c1dag, V_c2dag, V_c3, V_c4, kinetic)
        dim = length(hilbert_space_array)
        row_indices, col_indices, data = hamiltonian_sparse_data(data_length, hilbert_space_array, hilbert_space_dict, ts, t_c1dag, t_c2, Vs, V_c1dag, V_c2dag, V_c3, V_c4, kinetic)
        H = sparse(row_indices, col_indices, data, dim, dim)
        H
    end

end