for i = 1:nvb
   for j = 1:nent
     npost1 = vb_index(i); % index of post-synaptic neuron 1
     rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
     rpre1 = rf(C(npost1,:));
     [rpre2, rpost2] = unvec(j, 50);
     for m = 1:length(rpre1)
         dist = (rpost1 - rpost2)^2+(rpre1(m) - rpre2)^2;
         hx(i,j) = hx(i,j) + sigma^2*exp(-l*dist/2)*rpre1(m);
     end
   end
end
%upper left
hh = K+diag(seps.^2);

%lower left
xh = hx';

K = hh;
Ks = hx;
% Kss = xx;


% calculate posterior
K = preprocess(K, false);
L = chol(nearestSPD(K))';
alpha = L'\(L\(h - noise_mu));
%     v = L\Ks;
%     K_pos = Kss - v'*v; %posterior variance 
mu_pos =Ks'*alpha;
ret = reshape(mu_pos, nf, nf);


tmp = WRec_Novel*Diff_Rate;
mean(tmp)
mean(tmp(vb_index))
