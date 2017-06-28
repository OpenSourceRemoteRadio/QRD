 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  
 % Copyright (c) 2017, BigCat Wireless Pvt Ltd
 % All rights reserved.
 % 
 % Redistribution and use in source and binary forms, with or without
 % modification, are permitted provided that the following conditions are met:
 % 
 %     * Redistributions of source code must retain the above copyright notice,
 %       this list of conditions and the following disclaimer.
 %
 %     * Redistributions in binary form must reproduce the above copyright
 %       notice, this list of conditions and the following disclaimer in the
 %       documentation and/or other materials provided with the distribution.
 %
 %     * Neither the name of the copyright holder nor the names of its contributors
 %       may be used to endorse or promote products derived from this software
 %       without specific prior written permission.
 % 
 % 
 % 
 % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

matrix_size                                             = 16;

% For vectorization 4
vectorization                                           = 2;
qrd_instance                                            = 4;

input_matrix                                            = single(complex((rand(matrix_size,matrix_size)-1/2)*2,(rand(matrix_size,matrix_size)-1/2)*2));
input_matrix                                            = single(input_matrix);


%Column MAJOR
input_matrix_1d                                         = reshape((input_matrix),numel(input_matrix),1);
input_matrix_1d_uint32_real                             = typecast(real(single(input_matrix_1d)),'uint32');
input_matrix_1d_uint32_imag                             = typecast(imag(single(input_matrix_1d)),'uint32');
input_matrix_vectorized_real                            = transpose(reshape(input_matrix_1d_uint32_real,vectorization,numel(input_matrix_1d_uint32_real)/vectorization));
input_matrix_vectorized_imag                            = transpose(reshape(input_matrix_1d_uint32_imag,vectorization,numel(input_matrix_1d_uint32_imag)/vectorization));

input_matrix_vectorized_1_real                          = input_matrix_vectorized_real(:,1);
input_matrix_vectorized_1_imag                          = input_matrix_vectorized_imag(:,1);
input_matrix_vectorized_2_real                          = input_matrix_vectorized_real(:,2);
input_matrix_vectorized_2_imag                          = input_matrix_vectorized_imag(:,2);

input_matrix_vectorized_1_real                          = dec2hex(input_matrix_vectorized_1_real,8);
input_matrix_vectorized_1_imag                          = dec2hex(input_matrix_vectorized_1_imag,8);
input_matrix_vectorized_2_real                          = dec2hex(input_matrix_vectorized_2_real,8);
input_matrix_vectorized_2_imag                          = dec2hex(input_matrix_vectorized_2_imag,8);

fileID = fopen('../sim/input_data.txt','w+');

for index = 1:(matrix_size*matrix_size/vectorization)
        fprintf(fileID,'%s%s%s%s\n',input_matrix_vectorized_2_imag(index,:),input_matrix_vectorized_2_real(index,:),input_matrix_vectorized_1_imag(index,:),input_matrix_vectorized_1_real(index,:));
end
fclose(fileID);

n                                                       = matrix_size;
q                                                       = zeros(n,n);
r                                                       = zeros(n,n);
A_mod                                                   = input_matrix;
A_mod1                                                  = input_matrix;

index1                                                  = 1;
index2                                                  = 1;
index3                                                  = 1;
for col=1:n
    r2_acc                                              = 0;
    for row=1:n
        r2_acc                                          = r2_acc + (A_mod(row,col)'*A_mod(row,col));
    end
    r2(col,col)                                         = r2_acc;
    inv_r(col,col)                                      = 1/sqrt(r2(col,col));
    expected_output_matrix2(index2)                     = inv_r(col,col);
    index2                                              = index2 + 1;
    inv_r2(col,col)                                     = inv_r(col,col)*inv_r(col,col);
    for col1 = (col+1):n
        rn_acc                                          = 0;
        for row=1:n
            rn_acc                                      = rn_acc + (A_mod(row,col)'*A_mod(row,col1));
        end
        rn(col,col1)                                    = rn_acc;
        r(col,col1)                                     = rn(col,col1)*inv_r(col,col);
        expected_output_matrix3(index3)                 = r(col,col1);
        index3                                          = index3 + 1;
     
    end
    r(col,col)                                          = r2(col,col)*inv_r(col,col);    
    
    for col2 = (col+1):n
        rn_over_r2(col,col2)                            = rn(col,col2)*inv_r2(col,col);
    end
    
    for col3 = (col+1):n
        for row=1:n
            A_mod(row,col3)                             = A_mod(row,col3) - (rn_over_r2(col,col3)*A_mod(row,col));
        end
    end
    
    for row=1:n
        q(row,col)                                      = A_mod(row,col)*inv_r(col,col);
        expected_output_matrix1(index1)                 = q(row,col);
        index1                                          = index1 + 1;
    end
end

[Q,R]                                                   = qr(input_matrix);
inv_input_matrix                                        = inv(input_matrix);

q_t                                                     = q';

inv_mat                                                 = zeros(n,n);
for col=1:n
     inv_mat(:,col) = bs(q_t(:,col),r);
end


expected_output_matrix1_1d                              = reshape(expected_output_matrix1,numel(expected_output_matrix1),1);
expected_output_matrix2_1d                              = reshape(expected_output_matrix2,numel(expected_output_matrix2),1);
expected_output_matrix3_1d                              = reshape(expected_output_matrix3,numel(expected_output_matrix3),1);


keyboard

if(0)
fileID                                                  = fopen('../sim/output_data1.txt','r');
output_matrix_string                                    = fscanf(fileID,'%s');
fclose(fileID);

output_matrix_hex                                       = transpose(reshape(output_matrix_string,32,[]));
output_matrix_real_1                                    = output_matrix_hex(:,25:32);
output_matrix_imag_1                                    = output_matrix_hex(:,17:24);
output_matrix_real_2                                    = output_matrix_hex(:,9:16);
output_matrix_imag_2                                    = output_matrix_hex(:,1:8);

output_matrix_real_1                                    = typecast(uint32(hex2dec(output_matrix_real_1)),'single');
output_matrix_imag_1                                    = typecast(uint32(hex2dec(output_matrix_imag_1)),'single');
output_matrix_real_2                                    = typecast(uint32(hex2dec(output_matrix_real_2)),'single');
output_matrix_imag_2                                    = typecast(uint32(hex2dec(output_matrix_imag_2)),'single');

output_matrix_1                                         = complex(output_matrix_real_1,output_matrix_imag_1);
output_matrix_2                                         = complex(output_matrix_real_2,output_matrix_imag_2);

output_matrix                                           = [output_matrix_1 output_matrix_2];
output_matrix_1d                                        = reshape(transpose(output_matrix),[],1);


figure,plot(real(expected_output_matrix1_1d),'b--');
hold on
plot(real(output_matrix_1d),'b');
plot(imag(expected_output_matrix1_1d),'r--');
plot(imag(output_matrix_1d),'r');
legend('Real Expected','Real','Imag Expected','Imag');
max_absolute_error = max(abs(abs(expected_output_matrix1_1d)-abs(output_matrix_1d)));
disp('Maximum absolute error');
max_absolute_error

fileID                                                  = fopen('../sim/output_data2.txt','r');
output_matrix_string                                    = fscanf(fileID,'%s');
fclose(fileID);


output_matrix_hex                                       = transpose(reshape(output_matrix_string,16,[]));
output_matrix_real                                      = output_matrix_hex(:,9:16);
output_matrix_imag                                      = output_matrix_hex(:,1:8);

output_matrix_real                                      = typecast(uint32(hex2dec(output_matrix_real)),'single');
output_matrix_imag                                      = typecast(uint32(hex2dec(output_matrix_imag)),'single');

output_matrix                                           = complex(output_matrix_real,output_matrix_imag);
output_matrix_1d                                        = output_matrix;
output_matrix_1d                                        = reshape(transpose(output_matrix),[],1);


figure,plot(real(expected_output_matrix2_1d),'b--');
hold on
plot(real(output_matrix_1d),'b');
plot(imag(expected_output_matrix2_1d),'r--');
plot(imag(output_matrix_1d),'r');
legend('Real Expected','Real','Imag Expected','Imag');
max_absolute_error = max(abs(abs(expected_output_matrix2_1d)-abs(output_matrix_1d)));
disp('Maximum absolute error');
max_absolute_error

fileID                                                  = fopen('../sim/output_data3.txt','r');
output_matrix_string                                    = fscanf(fileID,'%s');
fclose(fileID);

output_matrix_hex                                       = transpose(reshape(output_matrix_string,16,[]));
output_matrix_real                                      = output_matrix_hex(:,9:16);
output_matrix_imag                                      = output_matrix_hex(:,1:8);

output_matrix_real                                      = typecast(uint32(hex2dec(output_matrix_real)),'single');
output_matrix_imag                                      = typecast(uint32(hex2dec(output_matrix_imag)),'single');

output_matrix                                           = complex(output_matrix_real,output_matrix_imag);
output_matrix_1d                                        = output_matrix;
output_matrix_1d                                        = reshape(transpose(output_matrix),[],1);


figure,plot(real(expected_output_matrix3_1d),'b--');
hold on
plot(real(output_matrix_1d),'b');
plot(imag(expected_output_matrix3_1d),'r--');
plot(imag(output_matrix_1d),'r');
legend('Real Expected','Real','Imag Expected','Imag');
max_absolute_error = max(abs(abs(expected_output_matrix3_1d)-abs(output_matrix_1d)));
disp('Maximum absolute error');
max_absolute_error

end


fileID                                                  = fopen('../sim/output_data.txt','r');
output_matrix_string                                    = fscanf(fileID,'%s');
fclose(fileID);

output_matrix_hex                                       = transpose(reshape(output_matrix_string,64,[]));
output_matrix_real_1                                    = output_matrix_hex(:,57:64);
output_matrix_imag_1                                    = output_matrix_hex(:,49:56);
output_matrix_real_2                                    = output_matrix_hex(:,41:48);
output_matrix_imag_2                                    = output_matrix_hex(:,33:40);
output_matrix_real_3                                    = output_matrix_hex(:,25:32);
output_matrix_imag_3                                    = output_matrix_hex(:,17:24);
output_matrix_real_4                                    = output_matrix_hex(:,9:16);
output_matrix_imag_4                                    = output_matrix_hex(:,1:8);


output_matrix_real_1                                    = typecast(uint32(hex2dec(output_matrix_real_1)),'single');
output_matrix_imag_1                                    = typecast(uint32(hex2dec(output_matrix_imag_1)),'single');
output_matrix_real_2                                    = typecast(uint32(hex2dec(output_matrix_real_2)),'single');
output_matrix_imag_2                                    = typecast(uint32(hex2dec(output_matrix_imag_2)),'single');
output_matrix_real_3                                    = typecast(uint32(hex2dec(output_matrix_real_3)),'single');
output_matrix_imag_3                                    = typecast(uint32(hex2dec(output_matrix_imag_3)),'single');
output_matrix_real_4                                    = typecast(uint32(hex2dec(output_matrix_real_4)),'single');
output_matrix_imag_4                                    = typecast(uint32(hex2dec(output_matrix_imag_4)),'single');

output_matrix_1                                         = complex(output_matrix_real_1,output_matrix_imag_1);
output_matrix_2                                         = complex(output_matrix_real_2,output_matrix_imag_2);
output_matrix_3                                         = complex(output_matrix_real_3,output_matrix_imag_3);
output_matrix_4                                         = complex(output_matrix_real_4,output_matrix_imag_4);



output_matrix                                           = [];
for iteration                                           = 1:(n/qrd_instance)
    iteration_result                                    = [output_matrix_1 output_matrix_2 output_matrix_3 output_matrix_4];
    iteration_result                                    = iteration_result(((iteration-1)*n)+1:((iteration-1)*n)+n,:);
    output_matrix                                       = [output_matrix iteration_result] ;
end

output_matrix_1d                                        = reshape((output_matrix),[],1);
expected_output_matrix_1d                               = reshape(inv_mat(:,:),[],1); 



figure,plot(real(expected_output_matrix_1d),'b--');
hold on
plot(real(output_matrix_1d),'b');
plot(imag(expected_output_matrix_1d),'r--');
plot(imag(output_matrix_1d),'r');
legend('Real Expected','Real','Imag Expected','Imag');
max_absolute_error                                      = max(abs(abs(expected_output_matrix_1d)-abs(output_matrix_1d)));
disp('Maximum absolute error');
max_absolute_error