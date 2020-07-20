#'Plot dfrr functional principal components
#'
#'Plot a \code{fpca.dfrr} object.
#'
#'@details This function plots the functional principal components,  contour plot and 3d surface of
#' the kernel function.
#'
#' If the package \code{\link[ggplot2]{ggplot2}} is installed, the contour plot of
#'  the kernel function is produced by setting the argument \code{plot.contour=TRUE}.
#'  Some graphical parameters of the contour plot can be modified by setting the (optional) argument
#'  \code{plot.contour.pars}.
#'
#'  If the package \code{\link[plotly]{plotly}} is installed, the 3d surface  of
#'  the kernel function is produced by setting the argument \code{plot.3dsurface=TRUE}.
#'  Some graphical parameters of the 3d surface can be modified by setting the (optional) argument
#'  \code{plot.3dsurface.pars}.
#'
#'@inheritParams fitted.dfrr
#'@param plot.contour a \code{boolean} indicating whether to print the contour plot of the kernel function.
#'It requires the package \code{\link[ggplot2]{ggplot2}} to be installed. Defaults to FALSE.
#'@param plot.eigen.functions a \code{boolean} indicating whether to print the principal components/eigen-functions. Defaults to TRUE.
#'@param plot.3dsurface a \code{boolean} indicating whether to print the 3d surface plot of the kernel function.
#'It requires the package \code{\link[plotly]{plotly}} to be installed. Defaults to FALSE.
#'@param  plot.contour.pars a named list of graphical parameters passed to the function \code{\link[ggplot2]{ggplot}}.
#'
#'@param  plot.3dsurface.pars a named list of graphical parameters passed to the function \code{\link[plotly]{plot_ly}}.
#'
#'@examples
#' set.seed(2000)
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' dfrr_fit<-dfrr(Y~X,yind=time)
#' fpcs<-fpca(dfrr_fit)
#' plot(fpcs,plot.eigen.functions=TRUE,plot.contour=TRUE,plot.3dsurface=TRUE)
#'
#'@method plot fpca.dfrr
#'
#'@export
plot.fpca.dfrr <-
function(dfrr_fit,plot.eigen.functions=TRUE,plot.contour=FALSE,plot.3dsurface=FALSE,
                          plot.contour.pars=list(breaks=NULL,minor_breaks = NULL,
                                                 n.breaks = NULL,
                                                 labels = NULL,
                                                 limits = NULL,
                                                 colors=NULL,
                                                 xlab=NULL,ylab=NULL,title=NULL),
                          plot.3dsurface.pars=list(xlab=NULL,ylab=NULL,zlab=NULL,
                                                    title=NULL,colors=NULL)
                          ){

  attr(dfrr_fit,"standardized") ->standardized
  attr(dfrr_fit,"dfrr_fit")     ->dfrr_fit

  if(plot.contour){

    if(require(ggplot2)){


      contour.pars.default<-list(breaks= waiver(),minor_breaks =  waiver(),n.breaks = NULL,labels = waiver(),limits = NULL,
                                 colors=c("#0D0887FF","#350498FF","#5402A3FF","#7000A8FF","#8B0AA5FF","#A31E9AFF","#B93289FF",
                                          "#CC4678FF","#DB5C68FF","#E97158FF","#F48849FF","#FBA139FF","#FEBC2AFF","#FADA24FF","#F0F921FF"),
                                 xlab="s",ylab="t",title=NULL)

      plot.contour.pars<-match.args.list(plot.contour.pars,contour.pars.default)


      mfrow<-par("mfrow")
      par(mfrow=c(1,1))

      time2<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out = 100)
      basis<-basis.dfrr(dfrr_fit)
      E2<-t(fda::eval.basis(time2,basis))
      if(standardized)
        sigma0<-dfrr_fit$sigma_theta_std
      else
        sigma0<-dfrr_fit$sigma_theta
      cvs2<-t(E2)%*%sigma0%*%E2

      x<-expand.grid(time2,time2)
      x<-data.frame(s=x[,1],t=x[,2],z=c(cvs2))
      v <- ggplot(x, aes(x=s, y=t, z = z))+
        geom_contour_filled(show.legend = FALSE)+
        scale_x_continuous(breaks = plot.contour.pars$breaks,
                           minor_breaks = plot.contour.pars$minor_breaks,
                           n.breaks = plot.contour.pars$n.breaks,
                           labels  = plot.contour.pars$labels ,
                           limits  = plot.contour.pars$limits ) +
        scale_y_continuous(breaks = plot.contour.pars$breaks,
                           minor_breaks = plot.contour.pars$minor_breaks,
                           n.breaks = plot.contour.pars$n.breaks,
                           labels  = plot.contour.pars$labels ,
                           limits  = plot.contour.pars$limits)+
        scale_fill_manual(values = plot.contour.pars$colors)+
        theme_minimal()+ labs(x = plot.contour.pars$xlab, y=plot.contour.pars$xlab)

      if(!is.null(plot.contour.pars$title))
        v<-v+labs(title=plot.contour.pars$title)

      plot(v)

      par(mfrow=mfrow)
    }else{
      warning("2D Coutour plot needs package 'ggplot2' to be installed")
    }

  }

  if(plot.eigen.functions){
    coef<-coef(dfrr_fit,return.principal.components = TRUE,standardized = standardized)
    plot(coef)
  }

  if(plot.3dsurface){
    if(require(plotly,quietly =TRUE)){

      surface3d.pars.default<-list(xlab="s",ylab="t",zlab="k(s,t)",title=NULL,
                                   colors=c("#0D0887FF","#350498FF","#5402A3FF","#7000A8FF","#8B0AA5FF","#A31E9AFF","#B93289FF",
                                            "#CC4678FF","#DB5C68FF","#E97158FF","#F48849FF","#FBA139FF","#FEBC2AFF","#FADA24FF","#F0F921FF"))
      plot.3dsurface.pars<-match.args.list(plot.3dsurface.pars,surface3d.pars.default)
      colorsLim<-seq(dfrr_fit$range[1],dfrr_fit$range[2],
                     length.out = length(plot.3dsurface.pars$colors))
      cols<-list()
      for(i in 1:length(colorsLim)){
        cols[[i]]<-c(colorsLim[i],plot.3dsurface.pars$colors[i])
      }

      time2<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out = 100)
      basis<-basis.dfrr(dfrr_fit)
      E2<-t(fda::eval.basis(time2,basis))
      if(standardized)
        sigma0<-dfrr_fit$sigma_theta_std
      else
        sigma0<-dfrr_fit$sigma_theta
      cvs2<-t(E2)%*%sigma0%*%E2

      m <- list(l = 0,  r = 0,  b = 0,  t = 0,  pad = 0)

      fig <- plot_ly(x = time2,y=time2,z=cvs2)%>%
        add_surface(colorscale =cols)%>%
        layout(margin=m,showlegend = TRUE,
               scene=list(xaxis=list(title=plot.3dsurface.pars$xlab),
                          yaxis=list(title=plot.3dsurface.pars$ylab),
                          zaxis=list(title=plot.3dsurface.pars$zlab)))
      print(fig)
    }else{
      warning("3D surface plot needs package 'plotly' to be installed")
    }
  }



}