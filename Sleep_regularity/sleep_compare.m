clear all
clc
%% regular sleep
sleep_reg = 1;
t = 24*365*50;
[tsleep_reg, C_reg,x] = sleep_irregular(t,sleep_reg);
C1_reg = C_reg(:,1);
C2_reg = C_reg(:,2);
C3_reg = C_reg(:,3);
C4_reg = C_reg(:,4);
C5_reg = C_reg(:,5);
C6_reg = C_reg(:,6);
C9_reg = C_reg(:,7);
C10_reg = C_reg(:,8);
%% irregular sleep
sleep_reg = 0;
% t = 24*365*10;
[tsleep_irreg, C_irreg,x] = sleep_irregular(t,sleep_reg);
C1_irreg = C_irreg(:,1);
C2_irreg = C_irreg(:,2);
C3_irreg = C_irreg(:,3);
C4_irreg = C_irreg(:,4);
C5_irreg = C_irreg(:,5);
C6_irreg = C_irreg(:,6);
C9_irreg = C_irreg(:,7);
C10_irreg = C_irreg(:,8);
%% plotting the results
figure(1)
plot(x, C1_reg,'r--',x,C1_irreg, 'r', x, C4_reg, 'b--',x, C4_irreg, 'b', 'Linewidth', 2);
legend('Ab40_b_r_a_i_n_p_a_r_e_n_c_h_y_m_a_-_r_e_g','Ab40_b_r_a_i_n_p_a_r_e_n_c_h_y_m_a_-_i_r_r_e_g', 'Ab42_b_r_a_i_n_p_a_r_e_n_c_h_y_m_a_-_r_e_g','Ab42_b_r_a_i_n_p_a_r_e_n_c_h_y_m_a_-_i_r_r_e_g')
xlabel('Hours'), ylabel('Ab #')
hold on
%% ab perivas
figure(2)
plot(x, C2_reg, 'g',x, C2_irreg, 'r', x, C5_reg, 'b',x, C5_irreg, 'm', 'Linewidth', 2);
legend('Ab40_p_e_r_i_v_a_s_-_r_e_g','Ab40_p_e_r_i_v_a_s_-_i_r_r_e_g', 'Ab42_p_e_r_i_v_a_s_-_r_e_g','Ab42_p_e_r_i_v_a_s_-_i_r_r_e_g')
xlabel('Hours'), ylabel('Ab #')
hold on
%% accumulation
figure(4)
plot(x, C9_reg, 'r--',x, C9, 'r', x, C10_reg, 'b--',x, C10, 'b', 'Linewidth', 2);
legend('Ab40_a_c_c_b_-_r_e_g','Ab40_a_c_c_b_-_c_a_f_f', 'Ab42_a_c_c_b_-_r_e_g','Ab42_a_c_c_b_-_c_a_f_f')
xlabel('Hours'), ylabel('Ab #')
hold on
figure(6)
plot(x, C3_reg, 'r--',x, C3, 'r', x, C6_reg, 'b--',x, C6, 'b', 'Linewidth', 2);
legend('Ab40_a_c_c_p_v_-_r_e_g','Ab40_a_c_c_p_v_-_c_a_f_f', 'Ab42_a_c_c_p_v_-_r_e_g','Ab42_a_c_c_p_v_-_c_a_f_f')
xlabel('Hours'), ylabel('Ab #')
%% time
figure(7);
plot(tsleep_irreg);