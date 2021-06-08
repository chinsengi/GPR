function [K, Kl, Ksn, Ksig] = Kgen(sigma, l, seps_neuron, vb_index, C, rf) %#codegen
   nvb = length(vb_index);
   K = zeros(nvb);
   Ksig = zeros(nvb);
   Kl = zeros(nvb);
   Ksn = zeros(nvb);
   for i = 1:nvb
       for j = i:nvb
         npost1 = vb_index(i); % index of post-synaptic neuron 1
         npost2 = vb_index(j); % index of post-synaptic neuron 2
         rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
         rpost2 = rf(npost2); % firing rate of post-synaptic neuron 2
         rpre1 = rf(C(npost1,:));
         rpre2 = rf(C(npost2, :));
         for m = 1:length(rpre1)
             for n = 1:length(rpre2)
                 dist = (rpost1 - rpost2)^2+(rpre1(m) - rpre2(n))^2;
                 multi_factor = rpre1(m)*rpre2(n);
                 if dist == 0
                     K(i,j) = K(i,j)+seps_neuron^2*multi_factor;
                     Ksn(i,j) = Ksn(i,j)+2*seps_neuron*multi_factor;
                 end
                 K(i,j) = K(i,j) + sigma^2*exp(-dist/(2*l^2))*multi_factor;
                 Ksig(i,j) = Ksig(i,j)+ 2*sigma*exp(-dist/(2*l^2))*multi_factor;
                 Kl(i,j) = Kl(i,j) + sigma^2*dist*exp(-dist/(2*l^2))*multi_factor/l^3;
             end
         end
       end
   end
   K = K + K' - diag(diag(K));
   Kl = Kl + Kl' - diag(diag(Kl));
   Ksn = Ksn + Ksn' - diag(diag(Ksn));
   Ksig = Ksig+Ksig' - diag(diag(Ksig));
end