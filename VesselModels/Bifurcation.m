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
        B1_A1 = 0;
        B2_A2 = 0;
        B3_A3 = 0;
        A1 = 0;
        B1 = 0;
        oms = 0;
        Yeff = 0;
    end
    
    methods

        function obj = Bifurcation(Rs, Ls, betas, rho, as, RWK1, RWK2, CWK)
                [obj.Rs, obj.Ls, obj.as] = deal(Rs, Ls, as);
                [obj.betas, obj.rho] = deal(betas, rho);
                [obj.RWK1, obj.RWK2, obj.CWK] = deal(RWK1, RWK2, CWK);
        end

        function obj = backpropagate(obj, inp)
            
            switch (obj.type)
                case 2
                     obj.oms = inp;
                     [obj.B1_A1,obj.B2_A2,obj.B3_A3,obj.Yeff] = type_II_bifurcation(...
                         obj.Rs(1),obj.Rs(2),obj.Rs(3),...
                         obj.Ls(1),obj.Ls(2),obj.Ls(2),...
                         obj.RWK1(1),obj.RWK2(1),obj.CWK(1),...
                         obj.RWK1(2),obj.RWK2(2),obj.CWK(2),...
                         obj.betas(1),obj.betas(2),obj.betas(3),...
                         obj.rho,obj.rho,obj.rho,...
                         obj.oms,obj.as(1),obj.as(2),obj.as(3));
                case 3
                    obj.oms = inp.oms;
                    obj.B3_A3 = inp.B1_A1;
                    [obj.B1_A1,obj.B2_A2,obj.Yeff] = type_III_bifurcation(...
                        inp.Yeff,obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,obj.oms,obj.as(1),obj.as(2));
                case 5
                    obj.oms = inp.oms;
                    obj.B3_A3 = inp.B1_A1;
                    [obj.A1,obj.B1,obj.B2_A2] = type_V_bifurcation(inp.Yeff,...
                        obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,obj.oms,obj.as(1),obj.as(2));
            end
        end

        function [Q2out, P2out, Q3out, P3out, A2, A3] = forwardpropagate(obj, P1outi)
            [Q3out, P3out, A3] = vesselforward(P1outi,obj.Ls(3),obj.Rs(3),obj.as(3),obj.oms,obj.rho,obj.betas(3),obj.B3_A3);
    
            [Q2out, P2out, A2] = vesselforward(P1outi,obj.Ls(2),obj.Rs(2),obj.as(2),obj.oms,obj.rho,obj.betas(2),obj.B2_A2);
        end

        function ves = vessel(obj, n)
            A = 0;
            t = 3;
            switch (n)
                case 1
                    WKP = [0, 0, 0];
                    if (obj.type == 5)
                        A = obj.A1;
                        B_A = obj.B1 ./ obj.A1;
                        t = 5;
                    else
                        B_A = obj.B1_A1;
                    end
                case 2
                    if (obj.type == 2)
                        WKP = [obj.RWK1(1), obj.RWK2(1), obj.CWK(1)];
                    else
                        WKP = [obj.RWK1, obj.RWK2, obj.CWK];
                    end
                    t = 2;
                    B_A = obj.B2_A2;
                case 3
                    if (obj.type == 2)
                        [obj.RWK1(2), obj.RWK2(2), obj.CWK(2)];
                        t = 2;
                    else
                        WKP = [0, 0, 0];
                    end
                    B_A = obj.B3_A3;
            end
            ves = Vessel(obj.Rs(n), obj.Ls(n), obj.as(n), obj.betas(n), obj.rho, WKP);
            ves.B1_A1 = B_A;
            ves.A1 = A;
            ves.type = t;
            ves.oms = obj.oms;
        end
    end
end

