classdef Bifurcation
    %BIFURCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = 3;
        Rs;
        Ls;
        betas;
        rho;
        as;
        RWK1;
        RWK2;
        CWK;
    end
    
    methods

        function obj = Bifurcation(Rs, Ls, betas, rho, as, RWK1, RWK2, CWK)
                [obj.Rs, obj.Ls, obj.as] = deal(Rs, Ls, as);
                [obj.betas, obj.rho] = deal(betas, rho);
                [obj.RWK1, obj.RWK2, obj.CWK] = deal(RWK1, RWK2, CWK);
        end

        function outs = backpropagate(obj, omega, Yeff_3)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes hereB3_A3
            switch (obj.type)
                case 2
                     [B1_A1,B2_A2,B3_A3,Yeff_1] = type_II_bifurcation(...
                         obj.Rs(1),obj.Rs(2),obj.Rs(3),...
                         obj.Ls(1),obj.Ls(2),obj.Ls(2),...
                         obj.RWK1(1),obj.RWK2(1),obj.CWK(1),...
                         obj.RWK1(2),obj.RWK2(2),obj.CWK(2),...
                         obj.betas(1),obj.betas(2),obj.betas(3),...
                         obj.rho,obj.rho,obj.rho,...
                         omega,obj.as(1),obj.as(2),obj.as(3));
                     outs = struct('B1_A1', B1_A1, 'B2_A2', B2_A2,...
                         'B3_A3', B3_A3, 'Yeff_1', Yeff_1);
                case 3
                    [B1_A1,B2_A2,Yeff_1] = type_III_bifurcation(...
                        Yeff_3,obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,omega,obj.as(1),obj.as(2));
                    outs = struct('B1_A1', B1_A1, 'B2_A2', B2_A2,...
                          'Yeff_1', Yeff_1);
                case 5
                    [A1,B1,B2_A2] = type_V_bifurcation(Yeff_3,...
                        obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,omega,obj.as(1),obj.as(2));
                    outs = struct('B1', B1, 'A1', A1, 'B2_A2', B2_A2);
            end
        end

        function [Q2out, P2out, A2, Q3out, P3out, A3] = forwardpropagate(obj, B2_A2, B3_A3, P1outi, omega)
            [Q3out, P3out, A3] = vessel(P1outi,obj.Ls(3),obj.Rs(3),obj.as(3),omega,obj.rho,obj.betas(3),B3_A3);
    
            [Q2out, P2out, A2] = vessel(P1outi,obj.Ls(2),obj.Rs(2),obj.as(2),omega,obj.rho,obj.betas(2),B2_A2);
        end
    end
end

