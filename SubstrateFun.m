function [DF,PM]=SubstrateFun(Substrate,FabrySwitch)

switch Substrate
    case 'NiFeLayer'
        DF=[1,1];
        PM=[1,1];
    case 'Si'
        DF=[1,3.41^2];
        PM=[1,1];
%     case 'SiO2'
%         DF=[1,4];
%         PM=[1,1];
    case 'SiO2-Si'
        DF=[1,4];
        PM=[1,1];
    case 'SiLorenzo'
        DF=[1,3.45^2];
        PM=[1,1];
    case 'SiReverse'
        DF=[3.41^2,1,1];
        PM=[1,1,1];
    case 'SiO2Au'
        DF=[1,2.1^2];
        PM=[1,1];
    case 'SiO2'
        DF=[1,2.1^2];
        PM=[1,1];
    case 'SiO2Reverse'
        DF=[2.1^2,1,1];
        PM=[1,1,1];
    case 'Si-Au'
        DF=[1,3.41^2];
        PM=[1,1];
    case 'SiO2-Au'
        DF=[1,2^2];
        PM=[1,1];
    otherwise
        fprintf('Warning! Unknown substrate!\n');
        bla
end

switch FabrySwitch
    case 'on'
        switch Substrate
            case {'Air','Si','SiO2','SiLorenzo'}
                DF=[DF,1];
                PM=[PM,1];
            case 'SiO2-Si'
                DF=[DF,3.41^2];
                PM=[PM,1];
            case 'SiO2Au'
                DF=[DF,(630.137+640.947i)^2];
                PM=[PM,1];
            case 'Si-Au'
                DF=[DF,(630.137+640.947i)^2];
                PM=[PM,1];
            case 'SiO2-Au'
                DF=[DF,(630.137+640.947i)^2];
                PM=[PM,1];
        end
        
end