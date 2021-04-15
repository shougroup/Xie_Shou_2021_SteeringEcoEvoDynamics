function score = score_wells(product, bio_T, coeff1, coeff2,mode)
if mode == 0
    pdtz = zscore(product);
    bioz = zscore(bio_T);
    b = [coeff1; coeff2];
    score = -[pdtz, bioz] * b;
end
if mode == 1
    score = product;
end
if mode == 2
    score = product ./ bio_T;
end
if ~or(mode ==0, or(mode == 1, mode ==2))
    error('mode must be 0, 1, or 2')
end
end
