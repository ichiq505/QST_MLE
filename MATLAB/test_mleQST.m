clc
disp("*** mleQST 테스트를 시작합니다. ***")

num_qubits = 1;
shots = 1000;
X = get_randstate(num_qubits);
% disp("▼▼ 임의 생성된 양자상태 X로 측정 시뮬레이션을 합니다. shot 수 = " + shots + " ▼▼")
% tic,[F,P] = do_paulimeas(X,shots); toc
% disp("▼▼ 측정 시뮬레이션을 통해 얻은 확률데이터로 QST를 수행합니다. ▼▼")
% tic,[Y, Ydiffs] = mleQST(num_qubits,F,1e-4,50,[.3 .7]); elapsedY = toc
% disp('Frobenius norm   /   Trace distance')
% disp([norm(X-Y,'fro'), get_statediff(X,Y,'trace')])
% disp("▼▼ Born rule으로 계산한 이론적인 확률데이터로 QST를 수행합니다. ▼▼")
% tic,[Z, Zdiffs] = mleQST(num_qubits,P,1e-4,50,[.3 .7]); elapsedZ = toc
% disp('Frobenius norm   /   Trace distance')
% disp([norm(X-Z,'fro'), get_statediff(X,Z,'trace')])

disp("▼▼ 임의 생성된 양자상태 X로 측정 시뮬레이션을 합니다. shot 수 = " + shots + ". 텐서곱에 내장함수 kron을 사용합니다. ▼▼")
tic,[F,P] = do_paulimeas_kron(X,shots); meas_simul_T = toc; disp(['경과 시간은 ',num2str(meas_simul_T),'초입니다.'])
disp("▼▼ 측정 시뮬레이션을 통해 얻은 확률데이터로 QST를 수행합니다. ▼▼")
tic,[Yk, Ydiffsk] = mleQST_kron(num_qubits,F,1e-4,50,[.3 .7]); elapsedYk = toc; disp(['경과 시간은 ',num2str(elapsedYk),'초입니다.'])
disp('Frobenius norm   /   Trace distance   /   Fidelity')
disp([norm(X-Yk,'fro'), get_statediff(X,Yk,'trace'), Fidelity(X,Yk)])
disp("▼▼ Born rule으로 계산한 이론적인 확률데이터로 QST를 수행합니다. ▼▼")
tic,[Zk, Zdiffsk] = mleQST_kron(num_qubits,P,1e-4,50,[.3 .7]); elapsedZk = toc; disp(['경과 시간은 ',num2str(elapsedZk),'초입니다.'])
disp('Frobenius norm   /   Trace distance   /   Fidelity')
disp([norm(X-Zk,'fro'), get_statediff(X,Zk,'trace'), Fidelity(X,Zk)])


clf, hold on
% subplot(221), semilogy(Ydiffs,'.-'), grid on, title(['convergence vs. iteration (sampling): ',num2str(elapsedY),'seconds'])
% subplot(223), semilogy(Zdiffs,'.-'), grid on, title(['convergence vs. iteration (theoretical probabilities): ',num2str(elapsedZ),'seconds'])
subplot(211), semilogy(Ydiffsk,'.-'), grid on, title(['convergence vs. iteration (sampling) [kron]: ',num2str(elapsedYk),'seconds'])
subplot(212), semilogy(Zdiffsk,'.-'), grid on, title(['convergence vs. iteration (theoretical probabilities) [kron]: ',num2str(elapsedZk),'seconds'])

figure(2)
% subplot(221), 
subplot(221), bar3(real(X)), axis([0 3 0 3 -inf inf]), title('Re(\rho^{target})')
% xticks(1:16),xticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
% yticks(1:16),yticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
xticks(1:2),xticklabels({'0','1'});
yticks(1:2),yticklabels({'0','1'});
subplot(222), bar3(real(Yk)), axis([0 3 0 3 -inf inf]), title('Re(\rho^{tomography})')
% xticks(1:16),xticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
% yticks(1:16),yticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
xticks(1:2),xticklabels({'0','1'});
yticks(1:2),yticklabels({'0','1'});% subplot(223),
subplot(223),bar3(imag(X)), axis([0 3 0 3 -inf inf]), title('Im(\rho^{target})')
% xticks(1:16),xticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
% yticks(1:16),yticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
xticks(1:2),xticklabels({'0','1'});
yticks(1:2),yticklabels({'0','1'});
subplot(224), bar3(imag(Yk)), axis([0 3 0 3 -inf inf]), title('Im(\rho^{tomography})')
% xticks(1:16),xticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
% yticks(1:16),yticklabels({'0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'});
xticks(1:2),xticklabels({'0','1'});
yticks(1:2),yticklabels({'0','1'});