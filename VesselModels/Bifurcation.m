classdef Bifurcation
    %BIFURCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = 3; % 2, 3 or 5
        Rs; % Radii of the 3 vessels [m]
        Ls; % Lengths of the 3 vessels [m]
        betas; % Stiffness factors of the 3 vessels [m/N]
        as; % Taper angle of the 3 vessels [rad]
        rho; % Density of blood [kg/m^3]
        RWK1; % First Windkessel Resistance(s) [Pa * s / m^3]
        RWK2; % Second Windkessel Resistance(s) [Pa * s / m^3]
        CWK; % Windkessel Capacitance(s) [m^3 / Pa]
        B1_A1 = 0; % B1/A1
        B2_A2 = 0; % B2/A2
        B3_A3 = 0; % B3/A3
        A1 = 0;
        B1 = 0;
        Yeff = 0; % Effective admittance at the inlet of the parent vessel
    end
    
    methods

        function obj = Bifurcation(Rs, Ls, betas, rho, as, RWK1, RWK2, CWK)
                [obj.Rs, obj.Ls, obj.as] = deal(Rs, Ls, as);
                [obj.betas, obj.rho] = deal(betas, rho);
                [obj.RWK1, obj.RWK2, obj.CWK] = deal(RWK1, RWK2, CWK);
        end

        function obj = backpropagate(obj, oms, prev)
            % Backpropagation using omegas, oms and successive vessel
            % object, prev to determine B1/A1 and possibly B2/A2 and B3/A3
            switch (obj.type)
                case 2
                     [obj.B1_A1,obj.B2_A2,obj.B3_A3,obj.Yeff] = type_II_bifurcation(...
                         obj.Rs(1),obj.Rs(2),obj.Rs(3),...
                         obj.Ls(1),obj.Ls(2),obj.Ls(2),...
                         obj.RWK1(1),obj.RWK2(1),obj.CWK(1),...
                         obj.RWK1(2),obj.RWK2(2),obj.CWK(2),...
                         obj.betas(1),obj.betas(2),obj.betas(3),...
                         obj.rho,obj.rho,obj.rho,...
                         oms,obj.as(1),obj.as(2),obj.as(3));
                case 3
                    obj.B3_A3 = prev.B1_A1;
                    [obj.B1_A1,obj.B2_A2,obj.Yeff] = type_III_bifurcation(...
                        prev.Yeff,obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,oms,obj.as(1),obj.as(2));
                case 5
                    obj.B3_A3 = prev.B1_A1;
                    [obj.A1,obj.B1,obj.B2_A2] = type_V_bifurcation(prev.Yeff,...
                        obj.Rs(1),obj.Rs(2),obj.Ls(1),obj.Ls(2),...
                        obj.RWK1,obj.RWK2,obj.CWK,...
                        obj.betas(1),obj.betas(2),...
                        obj.rho,obj.rho,oms,obj.as(1),obj.as(2));
            end
        end

        function [Q2out, P2out, Q3out, P3out, A2, A3] = forwardpropagate(obj, oms, P1outi)
            % Forwardpropagation using omegas, oms and length coordinate
            % in vessel, s as well as the pressure at the inlet of the
            % parent vessel, P0outi to determine pressure and flow rate
            ves3 = obj.vessel(3);
            ves2 = obj.vessel(2);
            [Q3out, P3out, A3] = ves3.forwardpropagate(oms, ves3.s(ves3.L), P1outi);
    
            [Q2out, P2out, A2] = ves2.forwardpropagate(oms, ves2.s(ves2.L), P1outi);
        end

        function ves = vessel(obj, n)
            % Returns a vessel object for any of the 3 vessels in the
            % bifurcation based on its number n. Use n=1 to get the parent
            % vessel, or n=2 or n=3 for the daughters.
            A = 0;
            t = 3;
            switch (n)
                case 1
                    WKP = [0, 0, 0];
                    if (obj.type == 5)
                        A = obj.A1;
                        B_A = obj.B1 ./ obj.A1;
                        t = 4;
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
                        WKP = [obj.RWK1(2), obj.RWK2(2), obj.CWK(2)];
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
        end
    end
end

