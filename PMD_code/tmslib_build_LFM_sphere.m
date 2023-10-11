function LFMm = tmslib_build_LFM_sphere(origin, coils, spos)
% TMSLIB_BUILD_LFM_SPHERE builds lead-field matrix for MEG forward problem
%
% Usage
%   LFMm = TMSLIB_BUILD_LFM_SPHERE(origin, coils, sources)
% Where
%   origin is the origin of the sphere
%   coils is a set of coils in the MEGBEM-format
%   sources is a N × 3 matrix of source positions
%
% Builds a lead-field matrix for MEG forward problem for a spherical head
% model centered around origin. Does not perform the sanity check of
% whether the coil is inside the 'head' (the sources outside the minimum
% radius for the coil). Uses the Sarvas formula (1987).
%
% Assumes that there is variables QPinds and QW, does not support QtoC!
%
% Author Lari Koponen
% Version 2016-02-18

    if ~iscell(spos)
        spos = {spos};
    end
    
    Nm = length(spos);
    N = size(coils.QP, 1);
    
    QP = coils.QP - ones(N, 1) * origin;
    QPN = sqrt(dots(QP, QP));
    QN = coils.QN;
    QW = coils.QW;
    QPinds = coils.QPinds;
    
    Q = eye(3);
    LFMm = cell(Nm, 1);
    for k = 1 : Nm ,
        s = spos{k} - ones(size(spos{k}, 1), 1) * origin;
        LFM_full = zeros(size(QP, 1), 3 * size(s, 1));
        for j = 1 : size(s, 1) ,
            r0 = s(j, :);
            for i = 1 : 3 ,
                av = QP - ones(N, 1) * r0;
                an = sqrt(dots(av, av));
                F = an .* (QPN .* (an + QPN) - dots(r0, QP));
                grad_F = ((an .* an ./ QPN + dots(av, QP) ./ an + 2 * (an + QPN)) * [1 1 1]) .* QP - ...
                    ((an + 2 * QPN + dots(av, QP) ./ an)) * r0;
                B = ((1e-7 ./ (F .* F)) * [1 1 1]) .* (F * crosses(Q(i, :), r0) - (triples(Q(i, :), r0, QP) * [1 1 1]) .* grad_F);
                LFM_full(:, 3 * j + i - 3) = QW .* dots(B, QN);
            end
        end
        LFM = zeros(size(QPinds, 1), 3 * size(s, 1)); 
        for i = 1 : size(QPinds, 1) ,
            LFM(i, :) = sum(LFM_full(QPinds(i, 1) : QPinds(i, 2), :), 1);
        end
        LFMm{k} = LFM;
    end
end

function c = dots(a, b)
    c = a(:, 1) .* b(:, 1) + a(:, 2) .* b(:, 2) + a(:, 3) .* b(:, 3);
end

function c = crosses(a, b)
    c = [a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2) ...
        a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3) ...
        a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1)];
end

function d = triples(a, b, c)
    d = (a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2)) .* c(:, 1) + ...
        (a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3)) .* c(:, 2) + ...
        (a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1)) .* c(:, 3);
end
