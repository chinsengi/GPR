function [K, KlE, KlI, KsnE, KsnI, KsigE, KsigI, gtg] = Kgen(param, vb_index, C, rf, rfpreI, rfpreE, postgrad) 
   nvb = numel(vb_index);
   nnE = numel(rfpreE);
   nnI = numel(rfpreI);
   nn = nnE+nnI;
   [sigmaE, lE, snE,...
              sigmaI, lI, snI, innoise] = extract_param(param);
   K = zeros(nvb);
   KsigE = zeros(nvb);
   KlE = zeros(nvb);
   KsnE = zeros(nvb);
   KsigI = zeros(nvb);
   KlI = zeros(nvb);
   KsnI = zeros(nvb);
   gtg = zeros(nvb);
   for i = 1:nvb
       for j = i:nvb
            npost1 = vb_index(i); % index of post-synaptic neuron 1
            npost2 = vb_index(j); % index of post-synaptic neuron 2
            rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
            rpost2 = rf(npost2); % firing rate of post-synaptic neuron 2
            rpre1I = rfpreI(C(npost1,nnE+1:nn));
            rpre1E = rfpreE(C(npost1,1:nnE));
            rpre2I = rfpreI(C(npost2,nnE+1:nn));
            rpre2E = rfpreE(C(npost2,1:nnE));
            for m = 1:length(rpre1I)
              for n = 1:length(rpre2I)
                 dist = (rpost1 - rpost2)^2+(rpre1I(m) - rpre2I(n))^2;
                 multi_factor = rpre1I(m)*rpre2I(n);
                 if dist == 0
                     K(i,j) = K(i,j)+snI^2*multi_factor;
                     KsnI(i,j) = KsnI(i,j)+2*snI*multi_factor;                        
                 end
                 K(i,j) = K(i,j) + sigmaI^2*exp(-dist/(2*lI^2))*multi_factor;
                 KsigI(i,j) = KsigI(i,j)+ 2*sigmaI*exp(-dist/(2*lI^2))*multi_factor;
                 KlI(i,j) = KlI(i,j) + sigmaI^2*dist*exp(-dist/(2*lI^2))*multi_factor/lI^3;
              end
            end
            for m = 1:length(rpre1E)
              for n = 1:length(rpre2E)
                 dist = (rpost1 - rpost2)^2+(rpre1E(m) - rpre2E(n))^2;
                 multi_factor = rpre1E(m)*rpre2E(n);
                 if dist == 0
                     K(i,j) = K(i,j)+snE^2*multi_factor;
                     KsnE(i,j) = KsnE(i,j)+2*snE*multi_factor;                       
                 end
                 K(i,j) = K(i,j) + sigmaE^2*exp(-dist/(2*lE^2))*multi_factor;
                 KsigE(i,j) = KsigE(i,j)+ 2*sigmaE*exp(-dist/(2*lE^2))*multi_factor;
                 KlE(i,j) = KlE(i,j) + sigmaE^2*dist*exp(-dist/(2*lE^2))*multi_factor/lE^3;
              end
            end
       end      
   end
   K = K + K' - diag(diag(K)); 
   KlE = KlE + KlE' - diag(diag(KlE)); 
   KsnE = KsnE + KsnE' - diag(diag(KsnE)); 
   KsigE = KsigE+KsigE' - diag(diag(KsigE)); 
   KlI = KlI + KlI' - diag(diag(KlI)); 
   KsnI = KsnI + KsnI' - diag(diag(KsnI)); 
   KsigI = KsigI+KsigI' - diag(diag(KsigI)); 
   gtg = gtg+gtg' - diag(diag(gtg)); 
end