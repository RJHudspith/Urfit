# F_pi

## Symmetric sources and sinks

This example input file illustrates how to compute the pion decay constant f_\pi from a CORR-type file.

     TrajXY = NULL,DATA/Meson_2pt_00ckn%zu.bin
     TrajStep = 52,8,116
     TrajStat = 1
     TrajFitr = 14.9,32.1
     TrajDims = 32,32,32,64
     TrajGsGk = 15,15,PLUS_PLUS
     TrajMom = 0,0,0
     TrajRW = NULL

     TrajXY = NULL,DATA/Meson_2pt_00ckn%zu.bin
     TrajStep = 52,8,116
     TrajStat = 1
     TrajFitr = 14.9,32.1
     TrajDims = 32,32,32,64
     TrajGsGk = 7,15,PLUS_MINUS
     TrajMom = 0,0,0
     TrajRW = NULL

     FitDef = PPAA
     Fit_NM = 1,1
     FitCorr = CORRELATED
     FitSims = 0,1,2
     FitTol = 1E-8
     FitMin = LM

     Guess_0 = 0.115
     Guess_1 = 100
     Guess_2 = 10

     Analysis = Correlator
     FileType = Corr_File

     Resample = BootStrap
     Nboots = 1000
     CovDiv = false
     CovBal = false
     CovEva = 1E-15

     Graph = UDBB.agr
     Graph_X = t/a
     Graph_Y = C(t)

The important things to note are the lines

    TrajGsGk = 15,15,PLUS_PLUS
    
    TrajGsGk = 7,15,PLUS_MINUS

That will fit a (sink,source) \gamma_5,\gamma_5 and \gamma_t\gamma_5,\gamma_5 correlator (in QDP gamma matrix conventions) and symmetrise in time with the folding PLUS_PLUS (cosh) and PLUS_MINUS (sinh)

We will be doing a PPAA fit, which does

    y_0 = p[1]*p[1]*( exp( -p[0]*x ) + exp( -p[0]*(Lt-x)) )
    
    y_1 = p[2]*p[1]*( exp( -p[0]*x ) - exp( -p[0]*(Lt-x)) )	

Which is useful when we have symmetric sources and sinks, e.g. Stochastic Wall sources or point sources or whatever.

The amplitudes for the Pseudoscalar P (p[1]) and the Temporal axial (p[2]) and the mass (p[0]) are all we need to compute the decay constant. The code will do this by the formula

    f_\pi/Z_A = \sqrt( 2*p[2]*p[2]/(V*p[0]) )

Where V will be determined from

    TrajDims = 32,32,32,64

i.e. 32x32x32

## Asymmetric sources and sinks

We can also get f_\pi from e.g. Coulomb gauge fixed wall sources RBC-style using mostly the same code but we need the four correlators P^W P^L, A^W A^L, P^W A^L and A^W P^L (probably 3 would suffice come to think of it), this is nice as they are all local at the sink so a smeared source would work too. Anyway we then call the

     PP_AA

Fit routine. Once we have the fit the decay contant is computed in the same fashion.

These parameters

     CovDiv = false
     CovBal = false
     CovEva = 1E-15

Should not be played with. In principle the code supports removing SVD "eigenvalues" (you would have to compile that in as by default it just does a Cholesky decomposition) and normalizing the correlation matrix (don't do this).