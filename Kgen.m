function [K, Kl, Ksn, Ksig, gtg] = Kgen(sigma, l, seps_neuron, vb_index, C, rf, innoise, postgrad) 
   nvb = length(vb_index);
   K = zeros(nvb);
   Ksig = zeros(nvb);
   Kl = zeros(nvb);
   Ksn = zeros(nvb);
   gtg = zeros(nvb);
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
                     if nargin > 6 % means there is input noise
                         % calculate the nearest pre-computed point with gradient
                         x = round(rpre1(m)); y = round(rpost1);
                         if x <= 0, x = 1; end
                         if y <= 0, y = 1; end
                         if x>50, x = 50; end
                         if y>50, y = 50; end
                         tmp_index = (y-1)*50+x;
                         tmp = postgrad(1, tmp_index)^2 + postgrad(2, tmp_index)^2;
                         gtg(i,j) = gtg(i,j) + tmp*multi_factor;
                         K(i,j) = K(i,j) + innoise^2*multi_factor*tmp;
                     end                         
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
   gtg = gtg+gtg' - diag(diag(gtg));
end