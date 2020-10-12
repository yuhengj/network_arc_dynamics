% ICD derivation, adapted from Scheinost et al., 2012

function norm_deg = icd(alpha,beta,tau)
norm_deg = exp(-alpha*tau^(beta));