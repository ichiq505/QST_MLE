function pauliPOVM = get_pauliPOVM_kron(num_qubits)
   pauli = {[1 1;1 1]/2,[1 -1;-1 1]/2,[1 -1i;1i 1]/2,[1 1i;-1i 1]/2,[1 0;0 0],[0 0;0 1]};
   
   neffects = 6^num_qubits;
   dim = 2^num_qubits;
    
   pauliPOVM = cell(1,neffects);
   effect_idx_num6 = zeros(1,num_qubits);
   for effect_idx=1:neffects
       effect_index_base6 = dec2base(effect_idx-1,6,num_qubits); for ii=1:num_qubits,   effect_idx_num6(ii) = str2double(effect_index_base6(ii)); end
        
       dim_sub = 2;
       effect = zeros(dim); 
       effect(1:2,1:2) = pauli{effect_idx_num6(1)+1};
       for ii=2:num_qubits
          dim_sub = bitshift(dim_sub,1);
          effect(1:dim_sub,1:dim_sub) = kron(effect(1:bitshift(dim_sub,-1),1:bitshift(dim_sub,-1)),pauli{effect_idx_num6(ii)+1});
       end
       
       pauliPOVM{effect_idx} = effect / (3^num_qubits);
   end
end