function M=MatrixMultiplier(varargin)

N=length(varargin);

M=varargin{1};
for uu=2:N
    tmp=M;
    tmp2=varargin{uu};
    for ii=1:2
        for jj=1:2
            M{ii,jj}=0;
            for kk=1:2
                M{ii,jj}=M{ii,jj}+tmp2{ii,kk}.*tmp{kk,jj};
            end
        end
    end
end