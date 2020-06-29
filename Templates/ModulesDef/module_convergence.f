      if(convergence.le.1) then
           call calculate_convergence(i, nprompt, nconv,
     >        convergence, annualglobalT,tsum, tsum_old, ddeltaT_old,
     >        annualglobalA, fcTOT, iceTOT, pressPtot, sigmaCRIT, 
     >        sigmaLIM, Tmax, f, tts, t2, exitFLAG)

      endif
