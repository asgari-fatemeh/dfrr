AMCEM <-
  function(
    Y,
    time,
    xData=NULL, #Nxq matrix of covariates
    T_E=500, #MaxNumber of E_steps
    T_G=80, #Number of samples to estimate integrals at each step
    T_G_cv=20,
    reltol=0.005, #b,sigma0,sigma_2
    J=NA,
    lambda_tol_sd=0.05,
    times_to_evaluate=seq(0,1,length.out = 100),
    seed=NA,
    method=c("REML","ML"),
    quiet=FALSE,
    ids=NULL,
    rangeval=c(0,1),
    normalize_beta=FALSE,
    basis=NULL
  ){
    method<-match.arg(method)
    flat_run<-1
    acceptRate_flatRun<-0.5
    max_iters_after_rejRate_one<-Inf
    imax_count<-0

    pen_order<-2
    bspline_order<-4

    AMCEM_converged<-FALSE

    if(!is.na(seed))
      set.seed(seed)

    T_G2=T_G

    library(fda)
    library(MASS)
    library(refund)
    library(tmvtnorm)

    N<-length(Y)

    if(is.null(ids))
      ids<-1:N


    if(is.null(xData)){
      q<-0
      estimate_b<-FALSE
    }else{
      q<-dim(xData)[2]
      estimate_b<-TRUE
    }


    M<-c()
    for(i in 1:N)
      M[i]<-length(time[[i]])

    if(is.na(J))
      J<-max(M)*1-4

    rel_tol_loglombdaSlope<-0.01
    lambda_tol_sd<-0.1
    sigma_tol_sd<-0.01
    cor_sigma<-0.3
    cor_lambda<-0.5


    t_E_tol<-10
    tol_b<-reltol*1.5
    tol_sigma0<-reltol



    rejRate_flatRun<-acceptRate_flatRun

    #A remedy for identifiablity problem when the problem is illposed
    lambdas_smoothing<-rep(0.00,1)


    burn.in.samples.main<-100
    burn.in.samples.light<-10
    t_G_multiplication<-1.02

    Rejrates<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)
    #set.seed(2000)


    time_start<-Sys.time()
    lambda_fixed<-FALSE
    sigma_fixed<-FALSE
    b_fixed<-FALSE
    sigma0_fixed<-FALSE


    fix_rejRate_at_one<-FALSE
    lambda_step<-sigma_step<-b_step<-sigma0_step<-NA

    sigma2Counter<-0
    if (is.null(basis))
      basis<-create.bspline.basis(rangeval = rangeval,nbasis = J,norder=bspline_order)
    else
      J<-basis$nbasis



    #basis<-create.fourier.basis(nbasis = J)

    inee<-fda::inprod(basis,basis,rng=rangeval)

    #################################

    time2<-seq(rangeval[1],rangeval[2],length.out = 100)
    E2<-t(eval.basis(time2,basis))

    P2es<-list()
    P2esRes<-list()

    for(i in 1:N)
      P2es[[i]]<-eval.penalty(basis,pen_order,rng=c(max(rangeval[1],min(time[[i]])-0.02),
                                                    min(rangeval[2],max(time[[i]])+0.02)))



    if(estimate_b)
      b0<-t(t(rep(0,J*q)))
    else{
      b0<-t(t(matrix(0,2,J)))
      B<-B2<-matrix(0,2,J)
    }


    sigma0<-diag(nrow = J)#cov(Y)
    sigma_2<-0.5

    P2<-eval.penalty(basis,pen_order,rng=rangeval)

    P22<-fda::eval.penalty(basis,pen_order,rng=rangeval)
    P00<-diag(nrow=J)

    #P2<-alpha*P00+(1-alpha)*P22

    E<-list()
    X<-list()
    for (i in 1:N) {
      E[[i]]<-t(eval.basis(time[[i]],basis))
      if(M[i]==1)
        E[[i]]<-t(t(c(E[[i]])))
      if(estimate_b)
        X[[i]]<-kronecker(t(xData[i,]),diag(nrow = J))
    }



    if(estimate_b){
      I<-diag(nrow = N)
      XX<-xData
      A<-I-XX%*%solve(t(XX)%*%XX)%*%t(XX)
      U<-eigen(A)$vector
      U<-U[,-seq(N-q+1,N)]
    }else{
      U<-diag(nrow=N)
    }



    #z0<-rep(0,M)
    Z_tr<-matrix(0,N,M)

    lambdass<-list()
    Bs<-list()
    sigma0s<-list()
    sigma_2s<-c()
    lambda_EM<-list()
    lambdam<-c()
    zplist<-list()
    zflist<-list()
    par(mfrow=c(2,2))
    count<-0

    alpha<-0.0
    P2<-alpha*P00+(1-alpha)*P22


    for (t_E in 1:T_E) {
      #alpha<-0.1+0.9^(t_E-1)
      if(!quiet)
        pb<-txtProgressBar(0,N,style = 3)
      kk<-0
      S2<-matrix(0,N,M)
      Zsum<-matrix(0,M,M)
      Z_tr<-matrix(0,N,M)

      burn_in=round(T_G*(1/3))
      burn_in<-100
      thin<-1

      mu1<-list()
      mu2<-list()
      sigma<-list()
      sigmai<-list()
      Z_tr<-list()

      lambdas<-rep(Inf,N)
      lambdas_tmp<-lapply(1:N, function(i){c()})

      if(lambda_fixed)
      {
        if(T_G<T_G2)
          T_G<-T_G2
        else
          T_G<-ceiling(t_G_multiplication*T_G)

      }

      for (i in 1:N) {

        if(M[i]==1)
          kttt<-t(E[[i]])%*%sigma0%*%E[[i]]
        else
          kttt<-diag(diag(t(E[[i]])%*%sigma0%*%E[[i]]))

        if(sigma_2<0.01)
          sigma_2<-0.01

        cv<-sigma_2*kttt

        sigmai[[i]]<-sigma0-sigma0%*%E[[i]]%*%
          solve(t(E[[i]])%*%sigma0%*%E[[i]]+cv)%*%
          t(E[[i]])%*%sigma0

        if(estimate_b){
          mu1[[i]]<-X[[i]]%*%b0-sigma0%*%E[[i]]%*%
            solve(t(E[[i]])%*%sigma0%*%E[[i]]+cv)%*%
            t(E[[i]])%*%X[[i]]%*%b0

          mu2[[i]]<-sigma0%*%E[[i]]%*%
            solve(t(E[[i]])%*%sigma0%*%E[[i]]+cv)
        }else{
          mu1[[i]]<-matrix(0,nrow=J)

          mu2[[i]]<-sigma0%*%E[[i]]%*%
            solve(t(E[[i]])%*%sigma0%*%E[[i]]+cv)
        }


        sigma[[i]]<-t(E[[i]])%*%sigma0%*%E[[i]]+cv
      }

      zsum<-matrix(0,N,J)
      zusum<-matrix(0,J,J)
      epsilon<-0

      B<-matrix(b0,nrow = q,byrow = TRUE)

      start.value<-NULL
      burn.in.samples<-100



      zslist<-lapply(1:T_G, function(j){matrix(0,N,J)})

      lambda_EM[[t_E]]<-list()




      lambda_tol<-lambda_thresh<-lambda_log_slope<-Inf
      lambda_log_slope_base<-1

      sigma_tol<-sigma_thresh<-Inf
      sigma_cor<-lambda_cor<-1
      if(t_E>t_E_tol+4){
        lambda_tol<-sd(lambdam[(t_E-t_E_tol):(t_E-1)])
        sigma_tol<-sd(sigma_2s[(t_E-t_E_tol):(t_E-1)])

        lambda_cor<-abs(cor(lambdam[(t_E-t_E_tol+1):(t_E-1)],
                            lambdam[(t_E-t_E_tol):(t_E-1-1)]))
        sigma_cor<-abs(cor(sigma_2s[(t_E-t_E_tol+1):(t_E-1)],
                           sigma_2s[(t_E-t_E_tol):(t_E-1-1)]))

        lambda_log_slope<-abs(log(lambdam)[t_E-1]-log(lambdam)[t_E-t_E_tol-1])/t_E_tol
        lambda_log_slope_base<-abs(log(lambdam)[t_E_tol+2]-log(lambdam)[2])/(t_E_tol)

        #lambda_log_slope<-abs(mean(log(lambdam)[(t_E-3):(t_E-1)])-mean(log(lambdam)[(t_E-t_E_tol-3):(t_E-t_E_tol-1)]))/t_E_tol

      }


      if(t_E>=2){
        lambda_thresh<-lambdam[t_E-1]*lambda_tol_sd
        sigma_thresh<-sigma_2s[t_E-1]*sigma_tol_sd
      }
      print(paste0("lambda_slope: ",lambda_log_slope," , ",lambda_log_slope_base))
      if(lambda_cor<cor_lambda | lambda_log_slope<rel_tol_loglombdaSlope | lambda_log_slope/lambda_log_slope_base<0.1){
        if(!lambda_fixed)
        {
          time_lambda<-Sys.time()
          lambda_step<-t_E
        }
        lambda_fixed<-TRUE
      }else{
        lambda_fixed<-FALSE
      }

      if(sigma_cor<cor_sigma)
      {
        if(!sigma_fixed){
          time_sigma<-Sys.time()
          sigma_step<-t_E-1
        }
        sigma_fixed<-TRUE
      }else{
        sigma_fixed<-FALSE
      }

      rejRate<-rejRate_flatRun
      if(t_E<=flat_run){
        rejRate<-rejRate_flatRun
      }else{

        #Rejrates<-c(0.8,0.9,1)
        Ns<-sample(1:N,N/4)
        if(N<=60)
          Ns<-1:N

        sigma_2_cv<-c()
        EE_cv<-c()
        sigma_theta_cv<-list()
        for(r in 1:length(Rejrates))
        {
          rejrate<-Rejrates[r]

          zplist2<-list()
          zflist2<-list()
          zslist2<-lapply(1:T_G_cv,function(jj){matrix(0,N,J)})
          for(i in Ns){

            start.value<-zplist[[i]][1,]
            burn.in.samples<-10

            ########### Cross validation #########
            if(estimate_b)
              vnu0<-t(E[[i]])%*%X[[i]]%*%b0
            else
              vnu0<-matrix(0,nrow = M[i])

            lb<-rep(-Inf,M[i])
            ub<-rep(Inf,M[i])
            lb[Y[[i]]==1]<-0
            ub[Y[[i]]==0]<-0

            kk5<-0
            repeat{
              if(M[i]==1)
                zprimes<-rtmvnorm(T_G_cv,mean=c(vnu0),sigma=sigma[[i]],algorithm = "rejection",
                                  lower=lb,upper=ub)
              else
                zprimes<-rtmvnorm(T_G_cv,mean=c(vnu0),sigma=sigma[[i]],algorithm = "gibbs",
                                  lower=lb,upper=ub,burn.in.samples=burn.in.samples,thin=10,
                                  start.value=start.value)

              if(M[i]==1)
                zprimes<-matrix(zprimes,ncol=M[i])

              if(any(is.na(zprimes)))
                kk5<-kk5+1
              else
                break
              if(kk5>20)
                stop("unable to generate zprimes")
            }

            if(t_E>3)
              lambdai<-quantile(lambda_EM[[t_E-1]][[i]],rejrate)
            else
              lambdai<-Inf

            P2e<-P2es[[i]]

            zs<-t(sapply(1:dim(zprimes)[1], function(i1){
              mui<-c(mu1[[i]]+mu2[[i]]%*%t(t(zprimes[i1,])))
              z1<-NA
              zKz1<-NA
              zKz0<-NA
              countt<-0
              maxReject<-10
              repeat{
                countt<-countt+1
                z<-MASS::mvrnorm(1,mui,sigmai[[i]])

                zKz3<-t(z)%*%P2e%*%t(t(z))
                if(is.na(zKz0))
                  zKz0<-zKz3

                zKz3<-zKz3

                if(zKz3<=lambdai)
                  return(c(z,NA,zKz0))

                if(is.na(zKz1))
                {z1<-z;zKz1<-zKz3}
                if(zKz3<zKz1)
                {z1<-z;zKz1<-zKz3}

                if(countt>=maxReject)
                  return(c(z1,NA,zKz0))

              }

            }))


            zs<-zs[,1:J]

            zplist2[[i]]<-zprimes
            zflist2[[i]]<-zs

            for(j in 1:T_G_cv)
              zslist2[[j]][i,]<-zs[j,]
            ######################################
          }
          if(estimate_b){
            I<-diag(nrow = length(Ns))
            XX2<-xData[Ns,]
            A2<-I-XX2%*%solve(t(XX2)%*%XX2)%*%t(XX2)
            U2<-eigen(A2)$vector
            U2<-U2[,-seq(length(Ns)-q+1,length(Ns))]
          }else{
            U2<-diag(nrow = length(Ns))
          }


          zusum<-matrix(0,J,J)
          for(j in 1:T_G_cv)
          {
            zss<-zslist2[[j]][Ns,]
            zusum<-zusum+t(t(U2)%*%zss)%*%(t(U2)%*%zss)/T_G_cv
          }
          sigma_t2<-zusum/(length(Ns)-q)
          sigma_theta_cv[[r]]<-sigma_t2

          epsilon<-0
          for(j in Ns){
            kttt<-diag(t(E[[j]])%*%sigma_t2%*%E[[j]])^-1#diag(nrow = M[j])
            if(M[j]==1)
              kttt<-t(kttt)
            else
              kttt<-diag(kttt)


            epsilon<-epsilon+sum(diag((zplist2[[j]]-zflist2[[j]]%*%E[[j]])%*%kttt%*%
                                        t(zplist2[[j]]-zflist2[[j]]%*%E[[j]])))/T_G_cv
          }
          sigma_2_cv[r]<-epsilon/sum(M[Ns])


          Xs<-XX[Ns,]

          Zs_cv<-list()
          EE1<-0
          Ezt<-matrix(0,J,length(Ns))
          innerMatrix<-kronecker(Xs%*%solve(t(Xs)%*%Xs)%*%t(Xs),sigma_theta_cv[[r]])

          for(j in 1:T_G_cv){
            Zs_cv[[j]]<-t(sapply(Ns, function(ii){zflist2[[ii]][j,]}))
            EE1<-EE1+t(c(t(Zs_cv[[j]])))%*%innerMatrix%*%t(t(c(t(Zs_cv[[j]]))))
            Ezt<-Ezt+t(Zs_cv[[j]])

          }
          EE1<-EE1/T_G_cv
          EE2<-t(c(Ezt))%*%innerMatrix%*%t(t(c(Ezt)))

          EE_cv[r]<-EE1-EE2

        }

        cvals<-c()

        Xs<-XX[Ns,]

        for(r in 1:length(Rejrates)){
          cvals[r]<- (N-q)/2*log(det(sigma_theta_cv[[r]]))+sum(M[Ns])/2*log(sigma_2_cv[r])
          #cvals[r]<- cvals[r]+1/2*log(det(kronecker(t(Xs)%*%Xs,solve(sigma_theta_cv[[r]]))))
          cvals[r]<- cvals[r]+1/2*EE_cv[r]
        }



        print(cvals)

        rejRate<-Rejrates[which.min(cvals)]
        #rejRate<-Rejrates[which.min(sigma_2_cv)]


      }

      if(!quiet){
        print(paste0("T_G: ",T_G))
        print(paste0("RejRate: ",rejRate))
      }


      for(i in 1:N){
        if(t_E>2){
          start.value<-zplist[[i]][1,]
          burn.in.samples<-burn.in.samples.light
        }else{
          start.value<-NULL
          burn.in.samples<-burn.in.samples.main
        }

        kk<-kk+1
        if(!quiet)
          setTxtProgressBar(pb,kk)

        if(estimate_b)
          vnu0<-t(E[[i]])%*%X[[i]]%*%b0
        else
          vnu0<-matrix(0,nrow = M[i])

        lb<-rep(-Inf,M[i])
        ub<-rep(Inf,M[i])
        lb[Y[[i]]==1]<-0
        ub[Y[[i]]==0]<-0

        kk5<-0
        repeat{
          if(M[i]==1)
            zprimes<-rtmvnorm(T_G,mean=c(vnu0),sigma=sigma[[i]],algorithm = "rejection",
                              lower=lb,upper=ub)
          else
            zprimes<-rtmvnorm(T_G,mean=c(vnu0),sigma=sigma[[i]],algorithm = "gibbs",
                              lower=lb,upper=ub,burn.in.samples=burn.in.samples,thin=10,
                              start.value=start.value)

          if(M[i]==1)
            zprimes<-matrix(c(zprimes),ncol=M[i])



          if(any(is.na(zprimes)))
            kk5<-kk5+1
          else
            break
          if(kk5>20)
            stop("unable to generate zprimes")
        }


        if(t_E>flat_run)
          lambdai<-quantile(lambda_EM[[t_E-1]][[i]],rejRate)
        else
          lambdai<-Inf

        P2e<-P2es[[i]]

        zs<-t(sapply(1:dim(zprimes)[1], function(i1){
          mui<-c(mu1[[i]]+mu2[[i]]%*%t(t(zprimes[i1,])))
          z1<-zKz1<-NA
          zKz0<-NA
          countt<-0
          maxReject<-20
          repeat{
            countt<-countt+1
            z<-MASS::mvrnorm(1,mui,sigmai[[i]])

            zKz3<-t(z)%*%P2e%*%t(t(z))

            if(is.na(zKz3))
              next
            if(is.na(zKz0))
              zKz0<-zKz3

            if(zKz3<=lambdai)
              return(c(z,NA,zKz3))

            if(is.na(zKz1))
            {z1<-z;zKz1<-zKz3}
            if(zKz3<zKz1)
            {z1<-z;zKz1<-zKz3}

            if(countt>=maxReject)
              return(c(z1,NA,zKz3))
          }

        }))

        lambda_EM[[t_E]][[i]]<-zs[,J+2]

        lambdas[i]<-mean(zs[,J+2])
        zs<-zs[,1:J]


        for(j in 1:T_G)
          zslist[[j]][i,]<-zs[j,]

        #epsilon<-epsilon+sum((zprimes-zs%*%inee1%*%E[[i]])^2)/T_G

        zplist[[i]]<-zprimes
        zflist[[i]]<-zs

        zsum[i,]<-colMeans(zs)

      }

      zusum<-matrix(0,J,J)
      zpsum<-matrix(0,J,J)
      zxsum<-matrix(0,J,q)
      for(j in 1:T_G)
      {
        zusum<-zusum+t(t(U)%*%zslist[[j]])%*%(t(U)%*%zslist[[j]])
        zpsum<-zpsum+t(zslist[[j]])%*%zslist[[j]]
        if(estimate_b)
          zxsum<-zxsum+t(zslist[[j]])%*%XX
      }
      zpsum<-zpsum/T_G
      zxsum<-zxsum/T_G
      zusum<-zusum/T_G

      if(!quiet)
        print(lambdas)
      lambdam[t_E]<-mean(lambdas)

      #sigma_2<-epsilon/sum(M)


      if(method=="REML"){

        Sigma_t<-(1/(N-q))*zusum

        if(estimate_b){
          S1<-matrix(0,q*J,q*J)
          S3<-matrix(0,q*J,1)
          for (i1 in 1:N) {
            S1<-t(X[[i1]])%*%solve(Sigma_t)%*%X[[i1]]+S1
            S3<-t(X[[i1]])%*%solve(Sigma_t)%*%t(t(zsum[i1,]))+S3
          }

          if(t_E<=length(lambdas_smoothing))
            lambda<-lambdas_smoothing[t_E]
          else
            lambda<-tail(lambdas_smoothing,1)

          P<-kronecker(diag(nrow=q),P2*lambda)
          b0<-solve(S1+2*P)%*%S3

        }

      }else{

        if(estimate_b){
          b0<-zxsum%*%solve(t(XX)%*%XX)
          Sigma_t<-(zpsum-b0%*%(t(XX)%*%XX)%*%t(b0))/N
        }else{
          Sigma_t<-zpsum/N
        }
      }

      sigma0<-Sigma_t

      resids<-list()
      epsilon<-0
      for(i in 1:N){
        kttt<-diag(t(E[[i]])%*%Sigma_t%*%E[[i]])^-1
        if(M[i]==1)
          kttt<-t(kttt)
        else
          kttt<-diag(kttt)

        # kttt<-diag(nrow = M[i])

        resids[[i]]<-colMeans(zplist[[i]]-zflist[[i]]%*%E[[i]])

        epsilon<-epsilon+sum(diag((zplist[[i]]-zflist[[i]]%*%E[[i]])%*%kttt%*%
                                    t(zplist[[i]]-zflist[[i]]%*%E[[i]])))/T_G


      }
      sigma_2<-epsilon/sum(M)






      kttt<-diag(Sigma_t)^-0.5



      Sigma_t<-t(E2)%*%sigma0%*%E2
      #sigma0<-solve(E2%*%t(E2))%*%E2%*%diag(diag(Sigma_t)^-0.5)%*%Sigma_t%*%diag(diag(Sigma_t)^-0.5)%*%t(E2)%*%solve(E2%*%t(E2))
      sigma02<-solve(E2%*%t(E2))%*%E2%*%diag(diag(Sigma_t)^-0.5)%*%Sigma_t%*%diag(diag(Sigma_t)^-0.5)%*%t(E2)%*%solve(E2%*%t(E2))

      if(estimate_b){
        B<-matrix(b0,nrow = q,byrow = TRUE)
        cc<-sqrt(t(B[1,])%*%inee%*%t(t(B[1,])))

        if(cc<0.2)
          cc<-1

        if(!normalize_beta)
          cc<-1

        B<-B/cc
        sigma0<-sigma0/(cc^2)
        # B<-B%*%E2%*%diag(diag(Sigma_t)^-0.5)%*%t(E2)%*%solve(E2%*%t(E2))
        B2<-B%*%E2%*%diag(diag(Sigma_t)^-0.5)%*%t(E2)%*%solve(E2%*%t(E2))
        b2<-c(B2)
      }



      lambdass[[t_E]]<-lambdas
      Bs[[t_E]]<-B
      sigma0s[[t_E]]<-sigma0
      sigma_2s[t_E]<-sigma_2

      regs<-B%*%E2
      regs_std<-B2%*%E2
      b0<-t(t(c(t(B))))



      ############
      ### EigenFunctions
      Omega<-inee
      eg<-eigen(inee)
      Omega2<-eg$vectors%*%diag(eg$values^0.5)%*%t(eg$vectors)
      Omega2n<-eg$vectors%*%diag(eg$values^-0.5)%*%t(eg$vectors)
      sigm<-Omega2%*%sigma02%*%Omega2

      eg<-eigen(sigm)
      nus<-eg$values
      psis<-Omega2n%*%eg$vectors
      psisE2<-t(psis)%*%E2


      if(!quiet){
        col1<-'black'

        if(sigma0_fixed)
          col1<-'red'
        par(mfrow=c(2,2))


        #################
        ## Regs
        if(t_E>1)
        {
          col1<-col2<-col3<-'black'

          if(lambda_fixed)
            col1<-'red'
          if(sigma_fixed)
            col2<-'red'
          if(b_fixed)
            col3<-'red'

          plot(time2,regs[1,],'l',col=col3,main = "Unstandardized Reg 1")

          if(q>1)
            plot(time2,regs[2,],'l',col=col3,main = "Unstandardized Reg 2")
          else
            plot.new()

          plot(1:t_E,sigma_2s[1:t_E],'l',col=col2)
          plot(1:t_E,log(lambdam[1:t_E]),'l',col=col1)




          col1<-'black'

          if(b_fixed)
            col1<-'red'

          plot(time2,regs_std[1,],'l',col=col1,main = "Standardized Reg 1")
          if(q>1)
            plot(time2,regs_std[2,],'l',col=col1,main = "Standardized Reg 2")
          else
            plot.new()
          plot(time2,psisE2[1,],'l',col=col1,main = "Standardized PC 1")
          plot(time2,psisE2[2,],'l',col=col1,main = "Standardized PC 2")

        }
      }


      eval1<-eigen(sigma02)$values

      if(t_E>5){
        #cheking for convergence
        deltak<-time2[2]-time2[1]
        beta_new<-B2%*%E2
        beta_old<-b_pre%*%E2
        dif_b0<-mean(sqrt(deltak*rowSums((beta_new-beta_old)^2))) #Remanian Integral

        dif_sigma<-sqrt(sum((eval1-sigma_pre)^2)) #Hilbert Schmidth norm

        b_thr<-tol_b*sqrt(mean(b2^2))
        s_thr<-tol_sigma0*sqrt(sum(eval1^2) )

        if(!b_fixed & dif_b0<b_thr){
          b_fixed<-TRUE
          time_b<-Sys.time()
          b_step<-t_E
        }

        if(!sigma0_fixed & dif_sigma< s_thr){
          sigma0_fixed<-TRUE
          time_sigma0<-Sys.time()
          sigma0_step<-t_E
        }



      }


      b_pre<-B2
      sigma_pre<-eval1


      if(imax_count>max_iters_after_rejRate_one |
         all(sigma_fixed,sigma0_fixed,lambda_fixed,b_fixed)){
        AMCEM_converged<-TRUE
        break
      }
    }


    if(AMCEM_converged==TRUE){
      # convergence_steps<-max(lambda_step,sigma_step,b_step,sigma0_step)
      # convergence_seconds<-max(round(difftime(time_lambda,time_start,units = "secs"),1),
      #                          round(difftime(time_sigma,time_start,units = "secs"),1),
      #                          round(difftime(time_b,time_start,units = "secs"),1),
      #                          round(difftime(time_sigma0,time_start,units = "secs"),1)
      convergence_steps<-max(lambda_step,sigma_step)
      convergence_seconds<-max(round(difftime(time_lambda,time_start,units = "secs"),1),
                               round(difftime(time_sigma,time_start,units = "secs"),1)
      )
    }
    else{
      convergence_steps<-NA
      convergence_seconds<-NA
    }


    ###### Preparing outputs
    Omega<-inee
    eg<-eigen(inee)
    Omega2<-eg$vectors%*%diag(eg$values^0.5)%*%t(eg$vectors)
    Omega2n<-eg$vectors%*%diag(eg$values^-0.5)%*%t(eg$vectors)
    sigm<-Omega2%*%sigma0%*%Omega2
    eg<-eigen(sigm)
    nus<-eg$values
    psis<-Omega2n%*%eg$vectors

    Omega<-inee
    eg<-eigen(inee)
    Omega2<-eg$vectors%*%diag(eg$values^0.5)%*%t(eg$vectors)
    Omega2n<-eg$vectors%*%diag(eg$values^-0.5)%*%t(eg$vectors)
    sigm<-Omega2%*%sigma02%*%Omega2
    eg<-eigen(sigm)
    nus2<-eg$values
    psis2<-Omega2n%*%eg$vectors

    Et<-t(eval.basis(times_to_evaluate,basis))

    Beta<-B%*%Et
    Beta_std<-B2%*%Et

    if(!estimate_b)
      B<-B2<-Beta<-Beta_std<-NULL

    fitted_coefs<-zsum
    if(!is.null(ids))
      rownames(fitted_coefs)<-ids


    zsum_std<-zsum%*%E2%*%diag(diag(Sigma_t)^-0.5)%*%t(E2)%*%solve(E2%*%t(E2))

    fitted_coefs_std<-zsum_std
    if(!is.null(ids))
      rownames(fitted_coefs_std)<-ids

    res<-list(B=B,sigma_theta=sigma0,sigma_theta_std=sigma02,sigma_2=sigma_2,zHat=zsum,
              nus=nus,Theta=t(psis),
              B_std=B2,nus_std=nus2,Theta_std=t(psis2),
              converged=AMCEM_converged,
              convergence_steps=convergence_steps,
              convergence_seconds=convergence_seconds,
              times_to_evaluate=times_to_evaluate,
              range=rangeval,basis=basis,
              Beta=Beta,
              Beta_std=Beta_std,
              psis=t(psis)%*%Et,
              psis_std=t(psis2)%*%Et,resids=resids,
              ids=ids,
              varnames=colnames(xData),
              fitted_coefs=fitted_coefs,fitted_coefs_std=fitted_coefs_std,
              data=list(Y=Y,time=time,xData=xData))

    return(res)
  }
