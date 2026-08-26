#' @title get_pretty_name
#'  
#' @description
#' 
#'    Gets pretty name for the generated simulation columns - delimited by spaces instead of periods, and properly capitalized.
#'
#' @param col the column name to prettify   
#' @param wrap Default=F, if T add a line break (backslash n) every 15 characters (at the next space).
#' @param non.rand.percent.grep what to replace .non.random.percent with in the pretty name 
#' 
#' @return the pretty name
#'
#' @export
get_pretty_name <- function(col, wrap=F, non.rand.percent.grep='Rejection Rate') {
  rep <- ''
  if(grepl('.non.random.percent',col)){
    rep <- paste0(' ', non.rand.percent.grep)
  }
  parts <- unlist(lapply(strsplit(gsub('\\.',' ', gsub('.non.random.percent',rep,col)), ' ')[[1]], function(x) { paste0(toupper(substring(x,1,1)),substring(x,2))}))
  if(wrap==F) {
    return(paste(parts, collapse = ' '))
  }
  
  currlen <- 0
  ret <- ''
  for(p in parts) {
    currlen <- currlen + nchar(p[1])
    ret <- paste0(ret,p)
    if(currlen>15){
      ret <- paste0(ret,'\n')
      currlen <- 0
    } else {
      ret <- paste0(ret,' ')
    }
  }
  return(ret)

}

#' @title get_special_theme
#'  
#' @description
#' 
#'    Gets a special ggplot theme with a square aspect ratio and larger fonts.
#'  
#' @param lockAspectRatio DEFAULT=TRUE. If TRUE, lock square aspect ratio.
#'   
#' @return the ggplot2::theme
#'
#' @export
get_special_theme <- function(lockAspectRatio=TRUE) {
  if(lockAspectRatio){
      ggplot2::theme(aspect.ratio = 1, 
                     axis.title = ggplot2::element_text(size=12),
                     axis.text = ggplot2::element_text(size=10), 
                     legend.title = ggplot2::element_text(size=12),
                     legend.text = ggplot2::element_text(size=10))
  } else {
    ggplot2::theme(
                   axis.title = ggplot2::element_text(size=12),
                   axis.text = ggplot2::element_text(size=10), 
                   legend.title = ggplot2::element_text(size=12),
                   legend.text = ggplot2::element_text(size=10))
  }
}


#' @title plot_simulation_data
#' 
#' @description Outputs various plots (via ggplot2) for each Clumping Method (CM) in the columns of aggregated.nrp.heatmap.data, with the Degrees of Clumpedness (dfr) vs. the Current Sample Size (CSS).
#' 
#' @param aggregated.nrp.heatmap.data Dataframe result from calculate_nonrandom_percent_on_heatmap_data
#' @param wanted_ns relevant for doLinePlots==TRUE and/or doFlippedLinePlots==TRUE. If non-empty, only these values of Current Sample Size (css) will be displayed in these plots.
#'                  If empty (default), specific values of CSS will be used, as the following percentiles: 15, 25, 50, 75, 85, 100. 
#' @param wanted_ds relevant for doLinePlots==TRUE and/or doFlippedLinePlots==TRUE. If non-empty, only these values of Degrees Of Clumping (dfr) will be displayed in these plots.
#'                  If empty (default), all levels will be displayed.
#' @param use.nrp Default=TRUE, if TRUE the Non-Random Percents (ratio of Simulation Repetitions which acheived a value over the 95th percentile of the Large Random-Like Population Distribution) 
#'                for each CM are the ones used for the plots. Otherwise the numeric value of the CM is what is used.
#' @param doHeatmaps Default=TRUE, output a heatmap for each CM. The Fill (yellow to red) is the NRP or the mean CM value (depending on use.nrp), the y-axis is the levels of dfr, and the x-axis is the levels of CSS.
#' @param doLinePlots Default=FALSE, output a dfr-Comparison line plot for each CM. The y-axis is the NRP or the mean CM value (depending on use.nrp), the the x-axis is the levels of CSS, and there is a line for each level of dfr.
#' @param doFlippedLinePlots Default=FALSE, output a CSS-Comparison line plot for each CM. The y-axis is the NRP or the mean CM value (depending on use.nrp), the the x-axis is the levels of dfr, and there is a line for each level of CSS in wanted_ns
#' @param linewidth.lineplot DEFAULT=1, linewidth for both regular and flipped lineplots.
#' @param lineplot.min.color Default='blue', for min of gradient scaling of the colors of the lines.
#' @param lineplot.max.color Default='red', for max of gradient scaling of the colors of the lines.
#' @param AOI Default=NULL, if not-null, adds a rectangle around the AOI on the heatmap. Order must be: 
#' 
#' * min.css
#' * min.dfr
#' * max.css
#' * max.dfr
#'            
#' Note that it is safeguarded for order, and also that minimi and maximi are trimmed to the max possible value.
#' @param linewidth.for.AOI Default=1, if AOI is not null, the width of the rectangle of the AOI.
#' @param fname.in.title  Default=TRUE, if to display the pretty name of the quantification method as the main title or appended to the proper axis.
#' @param grep.pretty.name.non.random.percent pretty replacement for non.random.percent in plots
#' 
#' @return list of the ggplot2 plots
#' @export
plot_simulation_data <- function(aggregated.nrp.heatmap.data, wanted_ns = c(), wanted_ds = c(),
                                 use.nrp=TRUE, doHeatmaps=TRUE, doLinePlots=FALSE, doFlippedLinePlots=FALSE, linewidth.lineplot = 1, lineplot.min.color='blue', lineplot.max.color='red',
                                 AOI=NULL, linewidth.for.AOI=1, fname.in.title=TRUE, grep.pretty.name.non.random.percent='Rejection Rate'){
  ret <- list()
  if(length(wanted_ns)==0){
    wanted_ns<-quantile(aggregated.nrp.heatmap.data$css,c(0.15,0.25,0.5,0.75,0.85,1))
  }
  for(col in colnames(aggregated.nrp.heatmap.data)){
    if(xor(!use.nrp, grepl('.non.random.percent',col)) && !col %in% c('css','dfr')){
      pretty.name <- ifelse(fname.in.title,grep.pretty.name.non.random.percent,get_pretty_name(col, non.rand.percent.grep = grep.pretty.name.non.random.percent))
      main.title=ifelse(fname.in.title,get_pretty_name(col, non.rand.percent.grep = ''),'')
      if(doHeatmaps){
        css.offset <- max(aggregated.nrp.heatmap.data$css)-max(aggregated.nrp.heatmap.data[aggregated.nrp.heatmap.data$css<max(aggregated.nrp.heatmap.data$css),]$css)
        pretty.name.heatmap <- get_pretty_name(col, wrap=T, non.rand.percent.grep = grep.pretty.name.non.random.percent)
        filt <- aggregated.nrp.heatmap.data
        plt <- ggplot2::ggplot(filt,ggplot2::aes(y=dfr, x=css)) +
          ggplot2::geom_tile(ggplot2::aes(fill=.data[[col]])) +
          ggplot2::scale_fill_gradient(low="yellow", high="red", limits = c(0,1)) +
          ggplot2::xlim(1,max(filt$css)+css.offset/2) +
          ggplot2::labs(x='Current Sample Size (CSS)', y = 'Deviation From Random (DFR)', fill=ifelse(fname.in.title, grep.pretty.name.non.random.percent, pretty.name.heatmap), title=main.title) + 
          ggplot2::theme_classic() + get_special_theme()
        
        if(!is.null(AOI)) {
          if(length(AOI)!=4) {
            stop(paste('AOI has',length(AOI),'variables - expected 4!'))
          }
          css.min <- max(AOI[1], min(aggregated.nrp.heatmap.data$css))
          dfr.min <- max(AOI[2], min(aggregated.nrp.heatmap.data$dfr))
          css.max <- min(AOI[3], max(aggregated.nrp.heatmap.data$css))
          dfr.max <- min(AOI[4], max(aggregated.nrp.heatmap.data$dfr))
          
          if(css.min>=css.max) {
            stop(paste('AOI is invalid - css minimum exceeds maximum (',css.min,'>=',css.max,')'))
          }
          if(dfr.min>=dfr.max){
            stop(paste('AOI is invalid - dfr minimum exceeds maximum (',dfr.min,'>=',dfr.max,')'))
          }
          
          css.offset <- max(aggregated.nrp.heatmap.data$css)-max(aggregated.nrp.heatmap.data[aggregated.nrp.heatmap.data$css<max(aggregated.nrp.heatmap.data$css),]$css)
          dfr.offset <- max(aggregated.nrp.heatmap.data$dfr)-max(aggregated.nrp.heatmap.data[aggregated.nrp.heatmap.data$dfr<max(aggregated.nrp.heatmap.data$dfr),]$dfr)
          
          max.x.aoi <- css.max + css.offset/2
          min.x.aoi <- css.min - css.offset/2
          max.y.aoi <- dfr.max + dfr.offset/2
          min.y.aoi <- dfr.min - dfr.offset/2
          
          plt <- plt + ggplot2::geom_path(data=data.frame(x=c(max.x.aoi,min.x.aoi,min.x.aoi,max.x.aoi,max.x.aoi),
                                                          y=c(min.y.aoi,min.y.aoi,max.y.aoi,max.y.aoi,min.y.aoi)),
                                          ggplot2::aes(x=x,y=y),
                                          linewidth=linewidth.for.AOI)
        }
        
        ret <- append(ret,list(plt))
      }
      if(doLinePlots){
        filt <- aggregated.nrp.heatmap.data
        if(length(wanted_ds)>0){
          filt <- filt[filt$dfr %in% wanted_ds,]
        }
        plt <- ggplot2::ggplot(filt, ggplot2::aes(y=.data[[col]], x=css, col=dfr, group=dfr)) +
          ggplot2::geom_line(linewidth=linewidth.lineplot) + ggplot2::scale_color_gradient(low=lineplot.min.color,high = lineplot.max.color) +
          ggplot2::labs(x='Current Sample Size (CSS)', y=pretty.name, color="Deviation From Random (DFR)", title=main.title) + 
          ggplot2::theme_classic() + get_special_theme()
        
        if(use.nrp) {
          plt <- plt + ggplot2::ylim(0, 1)
        }
        
        ret <- append(ret,list(plt))
      }
      if(doFlippedLinePlots){
        plt <- ggplot2::ggplot(aggregated.nrp.heatmap.data[aggregated.nrp.heatmap.data$css %in% wanted_ns,], 
                               ggplot2::aes(y=.data[[col]], x=dfr, col=css, group=css)) +
          ggplot2::geom_line(linewidth=linewidth.lineplot) + ggplot2::scale_color_gradient(low=lineplot.min.color,high = lineplot.max.color) +
          ggplot2::labs(x='Deviation From Random (DFR)', y=pretty.name, color="Current Sample Size (CSS)", title=main.title) +
          ggplot2::theme_classic() + get_special_theme()
        
        if(use.nrp) {
          plt <- plt + ggplot2::ylim(0, 1)
        }
        
        ret <- append(ret,list(plt))
      }
    }
  }
  return(ret)
}

#' @title heatmap_control_plots
#' 
#' @description Outputs plots (via ggplot2) for each Clumping Method (CM), either comparing values of heatmap.data 
#' (simulations resulting from simulate_population_remove_from_original_for_smaller_pops) to the 
#' Large Random-Like Population Distributions (LRPD, result of generate_large_randomlike_population_distribution), or
#' Only displaying the LRPD with the Minimum Non-Random Threshold (MNRT, the 95th percentile of the LRPD) for simplicity.
#' 
#' @param heatmap.data dataframe returned from simulate_population_remove_from_original_for_smaller_pops, with columns per CM.
#' @param lrpd dataframe returned from generate_large_randomlike_population_distribution, with columns per CM.
#' @param wanted_n wanted population size on which to display the plots. Default 0 dictates to display for the largest size (Initial Sample Size, ISS)
#' @param just_hist Default=FALSE. If FALSE, display for each CM a histogram of the LRPD values, with a dashed blue line for the MNRT (95th percentile), and vertical lines of mean CM value from heatmap.data for each dfr.
#'                  If TRUE, display from each CM a histogram of the LRPD values, with a dashed blue line for the MNRT (95th percentile), such that values above it are red and values below are green. 
#' @param with_mnrt Default=FALSE. If TRUE and just_hist = TRUE, mnrt will be drawn as a dashed blue vertical line, and histogram values above will be red.
#'                  else if FALSE and just_hist = TRUE, mnrt will not be drawn, and hist vals will be green only.
#'                  If just_hist = FALSE, do nothing.
#' 
#' @return list of the ggplot2 plots
#' @export
heatmap_control_plots <- function(heatmap.data, lrpd, wanted_n = 0, just_hist=F, with_mnrt=F, mnrt.quantile=0.95) {
  
  ret <- list()
  if(wanted_n==0){
    wanted_n <- max(heatmap.data$css)
  }
  
  for(col in colnames(heatmap.data)[3:length(colnames(heatmap.data))]){
    if(just_hist){
      
      datarange <- c(min(lrpd[col]), max(lrpd[col]))
      
      if(with_mnrt){
          q95 <- quantile(lrpd[col],c(mnrt.quantile), na.rm=TRUE)
          plt <- ggplot2::ggplot(lrpd[lrpd[col]<=q95,], ggplot2::aes(x=.data[[col]])) + 
            ggplot2::geom_histogram(binwidth = (datarange[2]-datarange[1])/35 , fill='darkgreen') + 
            ggplot2::geom_histogram(data = lrpd[lrpd[col]>=q95,], binwidth = (datarange[2]-datarange[1])/35 , fill='firebrick') + 
            ggplot2::geom_vline(xintercept = q95, col='blue', linewidth=2, linetype='dashed') +
            ggplot2::theme_classic()
          
          ret <- append(ret,list(plt))
      } else {
        plt <- ggplot2::ggplot(lrpd, ggplot2::aes(x=.data[[col]])) + 
          ggplot2::geom_histogram(binwidth = (datarange[2]-datarange[1])/35 , fill='darkgreen') + 
          ggplot2::theme_classic()
        ret <- append(ret,list(plt))
      }
    }
    else {
      pts <- c()
      ds <- c()
      
      for(d in unique(heatmap.data$dfr)){
        sub <- heatmap.data[heatmap.data$css==wanted_n & heatmap.data$dfr==d,c(col)]
        pts <- c(pts, mean(sub))
        ds <- c(ds,c(d))
      }
      
      vline.data <- data.frame(percentage.clumped=factor(ds), x=pts, towards.uniform=as.factor(ds<0))
      
      
      datarange <- c(min(min(lrpd[col]),vline.data$x), max(max(lrpd[col]),vline.data$x))
      
      plt <- ggplot2::ggplot(lrpd, ggplot2::aes(x=.data[[col]])) + 
        ggplot2::geom_histogram(binwidth = (datarange[2]-datarange[1])/20 , fill='darkgreen') + 
        ggplot2::labs(title=paste(col,'- Random Distribution with dashed CI95,\n compared to means of different clumping levels, n=', wanted_n)) +
        ggplot2::geom_vline(xintercept = quantile(lrpd[col],c(mnrt.quantile), na.rm=TRUE), col='blue', linewidth=2, linetype='dashed') +
        ggplot2::xlim(datarange[1], datarange[2]) +
        ggplot2::geom_vline(data=vline.data, mapping=ggplot2::aes(xintercept=x, color=percentage.clumped, linetype=towards.uniform), linewidth=1.3) + 
        ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      ret <- append(ret,list(plt))
    }
  }
  return(ret)
}


#' @title plots_of_example_data
#'
#' @description given example data in example.data.with.values with a column for each Clumping Method (CM), output comparative plots for each CM, taking into account the Currrent Sample Size (CSS) of each population to see if relevant for the CM.
#' 
#' @param example.data.with.values Dataframe outputted by calculate_values_example_data
#' @param lrpd Dataframe outputted by generate_large_randomlike_population_distribution
#' @param cutoffs a vector with the same length as the amount of CMs in the example data. This is so that only populations whose size
#'                permits usage of the CM (after checking to see at which population size this is not relevant with the plots of
#'                plot_simulation_data). so every population with n<cutoffs_i will NOT appear in the plot of CM_i=example.data.with.values_(i+2) (since the first two columns are POPULATION and CSS)
#' @param withSimData Default=FALSE. If TRUE, for each CM the values of populations relevant for the cutoff of each CM are displayed as vertical lines, overlaid on a histogram
#'                    of the LRPD values with the Minimum Non-Random Threshold (95th percentile of the LRPD values) marked in a dashed blue line.
#'                    If FALSE, for each CM the values of the population are displayed in a dot plot, where the X axis is simply sorted according to the CM value on the y axis. The dots are sized according to the CSS of the population, and a shorthand of the first letter of each word in POPULATION is placed next to each dot.
#' @param withGradient default=FALSE. If withSimData=TRUE, it is ignored. If withGradient=FALSE, the plot as described above is displayed per CM. 
#'                     If withGradient=TRUE, it is overlain on a yellow-to-red gradient, and the value of numericCutoffForGradient comes into play, such that populations
#'                     with a CM value > numericCutoffForGradient will be colored purple, and those under green.
#' @param numericCutoffForGradient default=-1.7, which is the relevant cutoff for clark.evans.naive on the supplied bird.location.data. If withSimData=TRUE or withGradient=FALSE, it is ignored.
#'                                 If withSimData=FALSE and withGradient=TRUE, the aforementioned coloration is done, and the value appears as a black horizontal line.
#' @param line.width.gradient default=2, for gradient display, width of the vertical lines on the gradient.
#' @param hjust.labels default=0.5, hjust for labels of population names on the points
#' @param vjust.labels default=2, vjust for labels of population names on the points 
#' @param grep.pretty.name.non.random.percent pretty replacement for non.random.percent in legends
#' 
#' @return list of the ggplot2 plots
#' @export
plots_of_example_data <- function(example.data.with.values, lrpd, cutoffs, 
                                    withSimData=F, withGradient=F, 
                                    numericCutoffsForGradient=c(-2, -1),
                                    line.width.gradient = 2,
                                    hjust.labels = 0.5, vjust.labels = 2,
                                    grep.pretty.name.non.random.percent = 'Rejection Rate') {
  
  ret <- list()
  for(i in 4:length(colnames(example.data.with.values))){
    
    col <- colnames(example.data.with.values)[i]
    
    vline.data <- example.data.with.values[example.data.with.values$css>=cutoffs[i-3],c('POPULATION','css',col)]
    
    vline.data <- vline.data[order(vline.data[[col]]),]
    
    vline.data$shorthand <- ''
    
    
    for(j in 1:length(vline.data$shorthand)){
      res <- ''
      pts <- strsplit(vline.data[j,]$POPULATION,' ')[[1]]
      for(k in 1:length(pts)){
        p<-pts[k]
        if(k==length(pts)){
          ch <- p
        } else {
          ch <- substr(p,1,1)
        }
        res <- paste0(res,ch)
      }
      vline.data[j,]$shorthand <- res
    }
    
    if(withSimData){
      datarange <- c(min(min(lrpd[col]),min(vline.data[col])), max(vline.data[col]))
      
      plt <- ggplot2::ggplot(lrpd, ggplot2::aes(x=.data[[col]])) + 
        ggplot2::geom_histogram(binwidth = (datarange[2]-datarange[1])/20 , fill='darkgreen') + 
        ggplot2::geom_vline(xintercept = quantile(lrpd[col],c(0.95), na.rm=TRUE), col='blue', size=2, linetype='dashed') +
        ggplot2::xlim(datarange[1], datarange[2]) + ggplot2::xlab(get_pretty_name(col, non.rand.percent.grep = grep.pretty.name.non.random.percent)) + ggplot2::labs(color = 'Population') +
        ggplot2::geom_vline(data=vline.data, mapping=ggplot2::aes(xintercept=.data[[col]], color=POPULATION), size=1.3) + 
        ggplot2::theme_classic() + get_special_theme()
      
      ret <- append(ret,list(plt))
    } else {
      
      if(withGradient){
      
          color_function <- colorRampPalette(c('red', 'yellow'),alpha = TRUE)
          g <- grid::rasterGrob(color_function(100), width=grid::unit(1,"npc"), height = grid::unit(1,"npc"), interpolate = TRUE) 
          
          col.palette <- c()
          
          for(i in 1:nrow(vline.data)){
            col.palette[vline.data[i,1]] <- ifelse(vline.data[i,col] > numericCutoffsForGradient[2],'black',ifelse(vline.data[i,col] > numericCutoffsForGradient[1],'darkblue','darkgreen'))
          }
      }
      
      len <- length(vline.data$POPULATION)
      
      
      plt <- ggplot2::ggplot(vline.data, ggplot2::aes(x=1:length(POPULATION), y=.data[[col]], color=POPULATION, label=shorthand, size = css))
      
      if(withGradient){
        plt <- plt + ggplot2::annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      }
      
      plt <- plt +
        ggplot2::geom_point() +
        ggplot2::geom_text(hjust=hjust.labels,vjust=vjust.labels,size = 2.8, show.legend=FALSE)
      
      if(withGradient){
        plt <- plt + ggplot2::scale_color_manual(values=col.palette) +
          ggplot2::geom_hline(yintercept = numericCutoffsForGradient[1], color='black', linewidth = line.width.gradient) +
          ggplot2::geom_hline(yintercept = numericCutoffsForGradient[2], color='black', linewidth = line.width.gradient)
      }
        plt <- plt + ggplot2::xlim(0, len+1) +
          ggplot2::labs(x='Population (Sorted)', y=get_pretty_name(col, non.rand.percent.grep = grep.pretty.name.non.random.percent), color="Population", size = 'Current Sample Size (CSS)') +
          ggplot2::theme_classic() + get_special_theme(lockAspectRatio = !withGradient)
      
        ret <- append(ret,list(plt))
    }
  }
  
  return(ret)
}

