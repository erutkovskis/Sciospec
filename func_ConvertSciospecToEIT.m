% convert Sciospec data to EIT data

function Veit=func_ConvertSciospecToEIT(V,NChannel,NSkip,fulldata)

% V(i,j) : voltage between i-th electrode and the ground electrode subject to the j-th injection current
% V(:,j) : j-th column of the matrix V is the measured voltages with respect to j-th injection current

temp=eye(NChannel);
temp_new=[temp(:,NChannel) temp(:,1:(NChannel-1))];
Mat_Conv=eye(NChannel)-temp_new;

Veit_temp=Mat_Conv*(V);
Veit=Veit_temp(:);

if ~fulldata
    rmv_indx=func_rmv_skip(NChannel,NSkip);
    Veit(rmv_indx)=[];
end

% rmv_indx=[     
%     16     1     2    17    18    19    34    35    36    51    52    53    68    69    70    85    86    87   102   103   104   119   120   121
%    136   137   138   153   154   155   170   171   172   187   188   189   204   205   206   221   222   223   238   239   240   255   256   241
%    ]';

function rmv_indx=func_rmv_skip(N_channel,N_skip)

temp_inx=[ N_channel-N_skip:N_channel 1:(N_channel-N_skip-1) ; 1:N_channel ; (1+N_skip+1):N_channel 1:(1+N_skip)]';


temp_inx2=[temp_inx+[0:(N_channel-1)]'*N_channel*ones(1,3)]';

rmv_indx=temp_inx2(:);