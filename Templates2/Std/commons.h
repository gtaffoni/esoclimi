      common /tempmatrix/ tempmat
      common /gasplusvapor_pars/ pressPtot,cp_Ptot,molwtPtot
      common /zonalcloudcover/ zonalfc 
      common /interp_TOA/ Matrix_TOA,T_TOA,p_TOA,z_TOA,as_TOA 

      common /transportpar/ DTbc,T1bc,Sbc,vTgrad 
      common /meanZenDist/ mzd 

      common/matrices/f, fo,  fcmat, olrmat, albmat, asrmat,
     >     phi, dryfluxmat, atmfluxmat, tts, habmat, cxlmat, 
     >     boilmat, RGmat


      common/output1/Tmin, Tmax,
     >     printfhab, printchab, hcxl, TotOLR, sigmaRG,
     >     sigmaBoil, DelT_NE, meanTNorthHem, fhabN,
     >     meanAlbedoNorth, TotOLRN,
     >     fhab, chab,
     >     fom, t2
      common/output2/
     >     annualglobalT, DelT_EP, nhab,
     >     annualglobalA, TotASR, exitFLAG,
     >     fcTOT, iceTOT, phiDEG, zonalT, Htime, zonalOLR,
     >     zonalASR, zonalALB, zonalATFX, zonalDRFX
