function [N_input,N_output,indx] = draw_one_cell(N_input,N_output)
csN = cumsum(N_input);
indx = find(csN>=randi(csN(end)),1);
N_input(indx) = N_input(indx) - 1;
N_output(indx) = N_output(indx) + 1;
end