function pauliPOVM = get_pauliPOVM(num_qubits)
   pauli = {[1 1;1 1]/2,[1 -1;-1 1]/2,[1 -1i;1i 1]/2,[1 1i;-1i 1]/2,[1 0;0 0],[0 0;0 1]};
   pauli_vec = cell(size(pauli));
   for idx = 1:length(pauli_vec)
      pauli_vec{idx} = pauli{idx}(:); 
   end
   clear pauli idx
    
   neffects = 6^num_qubits;
   dim = 2^num_qubits;
    
   pauliPOVM = cell(1,neffects);
   effect_idx_num6 = zeros(1,num_qubits);
   col_num2 = zeros(1,num_qubits);
   row_num2 = zeros(1,num_qubits);
   for effect_idx=1:neffects
       effect_index_base6 = dec2base(effect_idx-1,6,num_qubits); for ii=1:num_qubits,   effect_idx_num6(ii) = str2double(effect_index_base6(ii)); end
        
       effect = ones(dim);
       for col = 1:dim
          col_base2 = dec2bin(col-1,num_qubits); for ii=1:num_qubits,  col_num2(ii) = str2double(col_base2(ii)); end
          for row = 1:dim
             row_base2 = dec2bin(row-1,num_qubits); for ii=1:num_qubits,   row_num2(ii) = str2double(row_base2(ii)); end
             row_num2 = row_num2 + 2*col_num2;
             for ii = 1:num_qubits
                effect(row,col) = effect(row,col)*pauli_vec{effect_idx_num6(ii)+1}(row_num2(ii)+1);
             end
          end
       end
       pauliPOVM{effect_idx} = effect / (3^num_qubits);
   end
end