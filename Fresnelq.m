%%------------Fresnel coefficient tp,ts as a function of a lateral wavevector (qx,qy)------------------
%lambda - wavelength of the calculated light
%DF - array of the dielectric function/permitivity (complex number) of the material from the top
%[superstrate (air), substrate1 (NiFe), substrate2 silicon ...)
%PM - permeability (usually array of ones)
%d - thickness of layers between superstrate and substrate (in our case
%usually only NiFe). Size will be smaller by 2 than PM/DF
%SourceLayerIndex - index of the layer where the source (induced
%polarization) is located. Usually it will be index 2
%OutputLayerIndex - index of the layer where we want to have our output
%(usually 1)
%---Outputs---
%htp - 
function [htp,hts]=Fresnelq(lambda,DF,PM,d,SourceLayerIndex,OutputLayerIndex)

% if nargin<7


% DF=[1,1.5,-5+0.5i,1.5,1];
% PM=[1,1.5,-5+0.5i,1.5,1];
% d=[30e-9,30e-9,Inf];
% SourceLayerIndex=3;
% lambda=500e-9;

htp=@tpFun;
hts=@tsFun;

N=length(DF);
kn=2*pi/lambda*sqrt(DF.*PM);


    function tp=tpFun(q)
        

        knz=cell(N,1);
        for ii=1:N
            knz{ii}=(kn(ii)^2-q.^2).^(1/2);
        end

        
        rpij=cell(N-1,1);
        tpij=cell(N-1,1);
        rpji=cell(N-1,1);
        tpji=cell(N-1,1);
        for ii=1:N-1
            tmp=DF(ii+1)*knz{ii};
            tmp2=DF(ii)*knz{ii+1};
            tmp3=tmp+tmp2;
            rpij{ii}=(tmp-tmp2)./tmp3;
            tpij{ii}=2*tmp./tmp3*sqrt((PM(ii+1)*DF(ii))/(PM(ii)*DF(ii+1)));
            rpji{ii}=-rpij{ii};
            tpji{ii}=2*tmp2./tmp3*sqrt((PM(ii)*DF(ii+1))/(PM(ii+1)*DF(ii)));
        end
        
        Pj=cell(length(d),1);
        PjInv=cell(length(d),1);
        for ii=2:N-1
            Pj{ii-1}=cell(2,2);
            if isinf(d(ii-1))==1
                rpij{ii}=0;
                %                 rpji{ii-1}=0;
                
                Pj{ii-1}{1,1}=1;
                Pj{ii-1}{1,2}=0;
                Pj{ii-1}{2,1}=0;
                Pj{ii-1}{2,2}=1;
                
                PjInv{ii-1}=Pj{ii-1};
            else
                tmp=exp(1i*knz{ii}*d(ii-1));
                Pj{ii-1}{1,1}=tmp;
                Pj{ii-1}{1,2}=0;
                Pj{ii-1}{2,1}=0;
                Pj{ii-1}{2,2}=1./tmp;
                
                PjInv{ii-1}{1,1}=Pj{ii-1}{2,2};
                PjInv{ii-1}{1,2}=0;
                PjInv{ii-1}{2,1}=0;
                PjInv{ii-1}{2,2}=Pj{ii-1}{1,1};
            end
        end
        

        
        Mij=cell(N-1,1);
        Mji=cell(N-1,1);
        for ii=1:N-1
            tmp=(tpij{ii}.*tpji{ii}-rpij{ii}.*rpji{ii});
            
            Mji{ii}=cell(2,2);
            Mji{ii}{1,1}=1;
            Mji{ii}{1,2}=-rpij{ii};
            Mji{ii}{2,1}=rpji{ii};
            Mji{ii}{2,2}=tmp;
            
            Mij{ii}=cell(2,2);
            Mij{ii}{1,1}=tmp;
            Mij{ii}{1,2}=rpij{ii};
            Mij{ii}{2,1}=-rpji{ii};
            Mij{ii}{2,2}=1;
        end
        
                
        MUp=cell(2,2);
        MUp{1,1}=1;
        MUp{1,2}=0;
        MUp{2,1}=0;
        MUp{2,2}=1;
        FactorUp=1;
        for ii=1:SourceLayerIndex-1
            MUp=MatrixMultiplier(MUp,Mji{ii},PjInv{ii});
            FactorUp=FactorUp.*tpji{ii};
        end
        
        MDown=Mij{end};
        FactorDown=tpij{end};
        for ii=N-2:-1:SourceLayerIndex
            MDown=MatrixMultiplier(MDown,Pj{ii},Mij{ii});
            FactorDown=FactorDown.*tpij{ii};
        end
        
        tmp=MDown{2,2}.*MUp{1,1}-MDown{1,2}.*MUp{2,1};

        switch OutputLayerIndex
            case 1
                tp=cell(1,2);
                tp{1,1}=FactorUp.*MDown{2,2}./tmp;
                tp{1,2}=FactorUp.*MDown{1,2}./tmp;
            case N
                tp=cell(1,2);
                tp{1,1}=FactorDown.*MUp{2,1}./tmp;
                tp{1,2}=FactorDown.*MUp{1,1}./tmp;
            case SourceLayerIndex
                tp=cell(2,2);
                tp{1,1}=MUp{2,1}.*MDown{1,2}./tmp;
                tp{1,2}=MUp{1,1}.*MDown{1,2}./tmp;
                tp{2,1}=MUp{2,1}.*MDown{2,2}./tmp;
                tp{2,2}=MUp{2,1}.*MDown{1,2}./tmp;
            otherwise
                if OutputLayerIndex<SourceLayerIndex
                    L=cell(2,2);
                    L{1,1}=1;
                    L{1,2}=0;
                    L{2,1}=0;
                    L{2,2}=1;
                    for ii=1:OutputLayerIndex-1
                        L=MatrixMultiplier(L,Mji{ii},PjInv{ii});
                    end

                    Factor=1;
                    for ii=OutputLayerIndex:SourceLayerIndex-1
                        Factor=Factor.*tpji{ii};
                    end

                    tp=cell(2,2);
                    tp{1,1}=Factor.*L{1,1}.*MDown{2,2}./tmp;
                    tp{1,2}=Factor.*L{1,1}.*MDown{1,2}./tmp;
                    tp{2,1}=Factor.*L{2,1}.*MDown{2,2}./tmp;
                    tp{2,2}=Factor.*L{2,1}.*MDown{1,2}./tmp;

                elseif OutputLayerIndex>SourceLayerIndex

                    L=Mij{end};
                    for ii=N-2:-1:OutputLayerIndex
                        L=MatrixMultiplier(L,Pj{ii},Mij{ii});
                    end

                    Factor=1;
                    for ii=OutputLayerIndex-1:-1:SourceLayerIndex
                        Factor=Factor.*tpij{ii};
                    end
                    
                    tp=cell(2,2);
                    tp{1,1}=Factor.*L{1,2}.*MUp{2,1}./tmp;
                    tp{1,2}=Factor.*L{1,2}.*MUp{1,1}./tmp;
                    tp{2,1}=Factor.*L{2,2}.*MUp{2,1}./tmp;
                    tp{2,2}=Factor.*L{2,2}.*MUp{1,1}./tmp;
                end
        end

    end

    function ts=tsFun(q)
        
        knz=cell(N,1);
        for ii=1:N
            knz{ii}=(kn(ii)^2-q.^2).^(1/2);
        end
        
        rsij=cell(N-1,1);
        tsij=cell(N-1,1);
        rsji=cell(N-1,1);
        tsji=cell(N-1,1);
        for ii=1:N-1
            tmp=PM(ii+1)*knz{ii};
            tmp2=PM(ii)*knz{ii+1};
            
            tmp3=tmp+tmp2;
            
            rsij{ii}=(tmp-tmp2)./tmp3;
            tsij{ii}=2*tmp./tmp3;
            rsji{ii}=-rsij{ii};
            tsji{ii}=2*tmp2./tmp3;
        end
        
        Pj=cell(length(d),1);
        PjInv=cell(length(d),1);
        for ii=2:N-1
            Pj{ii-1}=cell(2,2);
            if isinf(d(ii-1))==1
                rsij{ii}=0;
                %                 rpji{ii-1}=0;
                
                Pj{ii-1}{1,1}=1;
                Pj{ii-1}{1,2}=0;
                Pj{ii-1}{2,1}=0;
                Pj{ii-1}{2,2}=1;
                
                PjInv{ii-1}=Pj{ii-1};
            else
                tmp=exp(1i*knz{ii}*d(ii-1));
                Pj{ii-1}{1,1}=tmp;
                Pj{ii-1}{1,2}=0;
                Pj{ii-1}{2,1}=0;
                Pj{ii-1}{2,2}=1./tmp;
                
                PjInv{ii-1}{1,1}=Pj{ii-1}{2,2};
                PjInv{ii-1}{1,2}=0;
                PjInv{ii-1}{2,1}=0;
                PjInv{ii-1}{2,2}=Pj{ii-1}{1,1};
            end
        end
        
        Mij=cell(N-1,1);
        Mji=cell(N-1,1);
        for ii=1:N-1
            tmp=(tsij{ii}.*tsji{ii}-rsij{ii}.*rsji{ii});
            
            Mji{ii}=cell(2,2);
            Mji{ii}{1,1}=1;
            Mji{ii}{1,2}=-rsij{ii};
            Mji{ii}{2,1}=rsji{ii};
            Mji{ii}{2,2}=tmp;
            
            Mij{ii}=cell(2,2);
            Mij{ii}{1,1}=tmp;
            Mij{ii}{1,2}=rsij{ii};
            Mij{ii}{2,1}=-rsji{ii};
            Mij{ii}{2,2}=1;
        end
        
        MUp=cell(2,2);
        MUp{1,1}=1;
        MUp{1,2}=0;
        MUp{2,1}=0;
        MUp{2,2}=1;
        FactorUp=1;
        for ii=1:SourceLayerIndex-1
            MUp=MatrixMultiplier(MUp,Mji{ii},PjInv{ii});
            FactorUp=FactorUp.*tsji{ii};
        end
        
        MDown=Mij{end};
        FactorDown=tsij{end};
        for ii=N-2:-1:SourceLayerIndex
            MDown=MatrixMultiplier(MDown,Pj{ii},Mij{ii});
            FactorDown=FactorDown.*tsij{ii};
        end
        
        tmp=MDown{2,2}.*MUp{1,1}-MDown{1,2}.*MUp{2,1};

        switch OutputLayerIndex
            case 1
                ts=cell(1,2);
                ts{1,1}=FactorUp.*MDown{2,2}./tmp;
                ts{1,2}=FactorUp.*MDown{1,2}./tmp;
            case N
                ts=cell(1,2);
                ts{1,1}=FactorDown.*MUp{2,1}./tmp;
                ts{1,2}=FactorDown.*MUp{1,1}./tmp;
            case SourceLayerIndex
                ts=cell(2,2);
                ts{1,1}=MUp{2,1}.*MDown{1,2}./tmp;
                ts{1,2}=MUp{1,1}.*MDown{1,2}./tmp;
                ts{2,1}=MUp{2,1}.*MDown{2,2}./tmp;
                ts{2,2}=MUp{2,1}.*MDown{1,2}./tmp;
            otherwise
                if OutputLayerIndex<SourceLayerIndex
                    L=cell(2,2);
                    L{1,1}=1;
                    L{1,2}=0;
                    L{2,1}=0;
                    L{2,2}=1;
                    for ii=1:OutputLayerIndex-1
                        L=MatrixMultiplier(L,Mji{ii},PjInv{ii});
                    end

                    Factor=1;
                    for ii=OutputLayerIndex:SourceLayerIndex-1
                        Factor=Factor.*tsji{ii};
                    end

                    ts=cell(2,2);
                    ts{1,1}=Factor.*L{1,1}.*MDown{2,2}./tmp;
                    ts{1,2}=Factor.*L{1,1}.*MDown{1,2}./tmp;
                    ts{2,1}=Factor.*L{2,1}.*MDown{2,2}./tmp;
                    ts{2,2}=Factor.*L{2,1}.*MDown{1,2}./tmp;

                elseif OutputLayerIndex>SourceLayerIndex

                    L=Mij{end};
                    for ii=N-2:-1:OutputLayerIndex
                        L=MatrixMultiplier(L,Pj{ii},Mij{ii});
                    end

                    Factor=1;
                    for ii=OutputLayerIndex-1:-1:SourceLayerIndex
                        Factor=Factor.*tsij{ii};
                    end

                    ts=cell(2,2);
                    ts{1,1}=Factor.*L{1,2}.*MDown{2,1}./tmp;
                    ts{1,2}=Factor.*L{1,2}.*MDown{1,1}./tmp;
                    ts{2,1}=Factor.*L{2,2}.*MDown{2,1}./tmp;
                    ts{2,2}=Factor.*L{2,2}.*MDown{1,1}./tmp;
                end
        end

    end
end