#' @title population_type
#'  
#' @description
#'        possible population type names - clumped, random, uniform, and subrandom.
#' @export
population_type <- list(
  clumped = 'clumped',
  random = 'random',
  uniform = 'uniform',
  subrandom = 'subrandom'
)

#' @title simulation_type
#'  
#' @description
#' Possible simulation type names:
#'      \itemize{
#'         \item{random.to.clumped}: {Population in the half-range of random to clumped (DFR 0 to 100)}
#'         \item{uniform.to.random}: {Population in the half-range of uniform to random (DFR -100 to 0)}
#'         \item{uniform.to.random.to.clumped}: {Population in the full range of uniform to random to clumped (DFR -100 to 100)}
#'        }
#' @export
simulation_type <- list(
  random.to.clumped = 'random to clumped',
  uniform.to.random = 'uniform to random',
  uniform.to.random.to.clumped = 'uniform to random to clumped'
)

#' @title uniform_random_generation_type
#'  
#' @description
#' Possible ways to generate the random-uniform mixed populations:
#'      \itemize{
#'         \item{generate.mix.in.full.polygon}: {The populations are generated in the same way as the random-clumped mixed populations - 
#'         CSS x abs(DFR) uniform points within the full polygon, and CSS x (1 - abs(DFR)) random points within the full polygon.}
#'         \item{split.polygon.and.generate.separately}: {The populations are generated in different subpolygons that make up the full polygon,
#'         by first splitting the polygon into two approximately equal sized polygons of DFR and (1-DFR) ratios of the full polygon.
#'         Then, the uniform subpopulation is generated within the DFR sub-polygon, and the random subpopulation within the (1-DFR) subpolygon.
#'         This generation in our opinion is more correct, as it allows the uniform sub-population to retain their distance, and not have the 
#'         random additions encroach the uniformity.}
#'        }
#' @export
uniform_random_generation_type <- list(
  generate.mix.in.full.polygon = 'generate mix in full polygon',
  split.polygon.and.generate.separately = 'split polygon and generate separately'
)





#' @title check_clumping_method_validity
#'  
#' @description
#' 
#'   Given list of sf polygons and a function name, check if it works as a clumping method (CM):
#'      \itemize{
#'         \item{}{Exists in memory}
#'         \item{}{Accepts parameters polys, dataframe, and n.for.calc}
#'         \item{}{Return a numeric value}
#'         \item{}{The numeric value is rising as clumpedness rises}
#'        }
#' 
#' @param polys List of sf polygons to check on the clumping methods
#' @param cm_to_check Clumping Method (CM) to check
#' @param sim.type simulation type parameter - to know which direction to check.
#' @return A list with a success parameter (boolean), and a message of the error if success == FALSE.
#' ``` {r}
#'      ## If method is a valid CM
#'      list(success = T, message = '')
#'      
#'      ## If method is NOT a valid CM
#'      list(success = F, message = 'AN ERROR OCCURRED')
#' ```
#'
#' @export
check_clumping_method_validity <- function(polys, cm_to_check, sim.type){
  func <- NULL
  
  ## check if function exists
  tryCatch({
    func <- get(cm_to_check)
  }, error = function(e) {
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,'does not exist in environment')))
  })
  
  ## check if null at this point sometimes get does not collapse
  if(is.null(func)){
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,'does not exist in environment')))
  }
  
  example.data <- data.frame(TAG=1:100,X=0,Y=0)
  
  clumped <- generate_clumped_population(100, polys, max_dist_from_centroid_divisor = 2)
  random <- generate_random_population(100, polys)
  uniform <- generate_uniform_population(100, polys)
  
  ## check if function accepts values as needed
  ret.clumped <- NULL
  tryCatch({
    ret.clumped <- func(clumped, polys, 0)
  }, error = function(e) {
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,
                                             'does not accept dataframe and list of polygons,',
                                             'or collapses', e)))
  })
  
  ## check if numeric
  if(!is.numeric(ret.clumped)){
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,
                                             'accepts dataframe and list of polygons,',
                                             'but returns non-numeric', ret.clumped)))
  }
  
  ret.uniform <- func(uniform, polys, 0)
  ret.random <- func(random, polys, 0)
  
  ## check if rising with clustering
  if(sim.type != simulation_type$uniform.to.random && ret.random > ret.clumped){
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,
                                             'should rise with clustering',
                                             'but clustered value', ret.clumped,
                                             '< random value', ret.random)))
  }
  
  if(sim.type != simulation_type$random.to.clumped && ret.random < ret.uniform){
    return(list(success = F, message = paste('Clumping Method (CM)',cm_to_check,
                                             'should drop with uniformity',
                                             'but uniform value', ret.uniform,
                                             '> random value', ret.random)))
  }
  
  return(list(success = T, message = ''))
  
  
  
}

#' @title generate_simulation_polygons
#'  
#' @description
#' 
#'   Given a dataframe of individual locations, calculate and return a list of length 1
#'   with the convexhull (chull, otherwise known as a minimum convex polygon or MCP) 
#'   of the points as an sf polygon. If use_bounding_rect == TRUE, return a bounding rectangle of the convex hull.
#' 
#' @param individual.locations.data.frame Dataframe of individual locations, expecting columns X and Y
#' @param use_bounding_rect Default TRUE=T, whether to return an chull/MCP of the locations, or a bounding rectangle of this mcp.
#' @return A list of sf polygons with length 1, which can either be the chull/MCP of individual.locations.data.frame if use_bounding_rect=FALSE,
#'         or a bounding rectangle of said chull/MCP if use_bounding_rect=TRUE.
#' 
#' @export
generate_simulation_polygons <- function(individual.locations.data.frame, use_bounding_rect = T) {
  ## chull returns the row indices
  mcp_inds <- chull(individual.locations.data.frame[,c('X','Y')])
  ## close the polygon by repeating the first point
  mcp_inds <- c(mcp_inds,mcp_inds[1])
  
  mcp <- individual.locations.data.frame[mcp_inds,c('X','Y')]
  
  sf <- sfheaders::sfc_polygon(mcp)
  
  ## if use_bounding_rect, turn into a bounding rectangle.
  if(use_bounding_rect){
    
    df.rect <- data.frame(X=rev(c(min(mcp[,1]), max(mcp[,1]), max(mcp[,1]), min(mcp[,1]),min(mcp[,1]))),
                          Y=rev(c(min(mcp[,2]), min(mcp[,2]), max(mcp[,2]), max(mcp[,2]),min(mcp[,2]))))
    
    sf <- sfheaders::sfc_polygon(df.rect)
  }
  
  return(list(sf))
}


#' @title split_poly
#'  
#' @description
#' 
#' Split polygon into 2 APPROXIMATELY equal sized polygons, s.t the first will be approximately ratio.first.poly of the bottom (if h>w) or left (if h<=w) of the polygon,
#' and the second will be the remaining (1-ratio.first.poly) of the top (if h>w) or right (if h<=w) of the polygon. 
#' In case of ratio.first.poly==0 the first polygon will be null and the second the full polygon,
#' and vice versa if ratio.first.poly==1. 
#' Note that currently the split is according to the median of the longer of the the two axes of the polygon,
#' and therefore not necessarily will the area be equal to ratio.first.poly or 1-ratio.first.poly. However, this is inconsequential to the calculations afterwards, since
#' the uniform subpopulation will still be highly uniform within its subpolygon.
#' 
#' @param poly Original polygon as an sf polygon
#' @param ratio.first.poly 0-1, the approximate percentage to give the first of the two polygons.
#' @return A list of five parameters:
#'      \itemize{
#'         \item{ratio.first.poly}: {same as input, for tracking purposes}
#'         \item{full.poly}: {same as input polygon, for tracking purposes}
#'         \item{poly.1}: {a sub-polygon of approximately ratio.first.poly of the large polygon, or NULL if ratio.first.polygon==0. Will be at the bottom if h>w in the original polygon,
#'          and at the left if h<=w in the original polygon}
#'         \item{poly.2}: {a sub-polygon of approximately 1-ratio.first.poly of the large polygon, or NULL if ratio.first.polygon==1. Will be at the top if h>w in the original polygon,
#'          and at the right if h<=w in the original polygon}
#'         \item{plot}: {plot of the algorithmic steps, if wanted.}
#'        }
#' 
#' @export
split_poly <- function(poly, ratio.first.poly) {
  coords <- sf::st_coordinates(poly)
  xrange <- c(min(coords[,1]),max(coords[,1]))
  yrange <- c(min(coords[,2]),max(coords[,2]))
  
  w <- xrange[2] - xrange[1]
  h <- yrange[2] - yrange[1]
  
  poly.1 <- NULL
  poly.2 <- NULL
  
  if(h>w) {
    cutoff <- ratio.first.poly * h
    line.data <- data.frame(x=xrange, y=yrange[1] + cutoff, linestring_id=1)
    line <- sfheaders::sf_line(line.data)
    line.intersect <- sf::st_coordinates(sf::st_intersection(poly,line))
    
    poly.1.coords <- coords[coords[,2]<=min(line.data$y),]
    poly.2.coords <- coords[coords[,2]>=min(line.data$y),]
    
    nrows <- ifelse(is.null(nrow(poly.1.coords)),1,nrow(poly.1.coords))
    
    if(nrows>0){
      poly.1.coords <- matrix(nrow=0, ncol=2)
      line.added <- F
      for(i in 1:nrow(coords)){
        if(coords[i,2]>min(line.data$y) && line.added==F){
          poly.1.coords <- rbind(poly.1.coords,line.intersect[,1:2])
          line.added <- T
        }
        if(coords[i,2]<=min(line.data$y)) {
          poly.1.coords <- rbind(poly.1.coords,coords[i,1:2])
        }
      }
      if(nrow(poly.1.coords)>=3 && ratio.first.poly>0){
        tryCatch({
          poly.1 <- generate_simulation_polygons(poly.1.coords, use_bounding_rect = F)[[1]]
        }, error = function(e) {
          warning(paste('Failed to create poly.1:',e))
        })
      }
    }
    nrows <- ifelse(is.null(nrow(poly.2.coords)),1,nrow(poly.2.coords))
    if(nrows>0){
      poly.2.coords <- matrix(nrow=0, ncol=2)
      line.added <- F
      for(i in 1:nrow(coords)){
        if(coords[i,2]<min(line.data$y) && line.added==F){
          poly.2.coords <- rbind(poly.2.coords,line.intersect[,1:2])
          line.added <- T
        }
        if(coords[i,2]>=min(line.data$y)) {
          poly.2.coords <- rbind(poly.2.coords,coords[i,1:2])
        }
      }
      if(nrow(poly.2.coords)>=3 && ratio.first.poly<100){
        tryCatch({
          poly.2 <- generate_simulation_polygons(poly.2.coords, use_bounding_rect = F)[[1]]
        }, error = function(e) {
          warning(paste('Failed to create poly.2:',e))
        })
      }
    }
    
  } else {
    cutoff <- ratio.first.poly * w
    line.data <- data.frame(x=xrange[1] + cutoff, y=yrange, linestring_id=1)
    line <- sfheaders::sf_line(line.data)
    line.intersect <- sf::st_coordinates(sf::st_intersection(poly,line))
    
    poly.1.coords <- coords[coords[,1]<=min(line.data$x),]
    poly.2.coords <- coords[coords[,1]>=min(line.data$x),]
    
    nrows <- ifelse(is.null(nrow(poly.1.coords)),1,nrow(poly.1.coords))
    
    if(nrows>0){
      poly.1.coords <- matrix(nrow=0, ncol=2)
      line.added <- F
      for(i in 1:nrow(coords)){
        if(coords[i,1]>min(line.data$x) && line.added==F){
          poly.1.coords <- rbind(poly.1.coords,line.intersect[,1:2])
          line.added <- T
        }
        if(coords[i,1]<=min(line.data$x)) {
          poly.1.coords <- rbind(poly.1.coords,coords[i,1:2])
        }
      }
      if(nrow(poly.1.coords)>=3 && ratio.first.poly>0){
        tryCatch({
          poly.1 <- generate_simulation_polygons(poly.1.coords, use_bounding_rect = F)[[1]]
        }, error = function(e) {
          warning(paste('Failed to create poly.1:',e))
        })
      }
    }
    nrows <- ifelse(is.null(nrow(poly.2.coords)),1,nrow(poly.2.coords))
    if(nrows>0){
      poly.2.coords <- matrix(nrow=0, ncol=2)
      line.added <- F
      for(i in 1:nrow(coords)){
        if(coords[i,1]<min(line.data$x) && line.added==F){
          poly.2.coords <- rbind(poly.2.coords,line.intersect[,1:2])
          line.added <- T
        }
        if(coords[i,1]>=min(line.data$x)) {
          poly.2.coords <- rbind(poly.2.coords,coords[i,1:2])
        }
      }
      if(nrow(poly.2.coords)>=3 && ratio.first.poly<100){
        tryCatch({
          poly.2 <- generate_simulation_polygons(poly.2.coords, use_bounding_rect = F)[[1]]
        }, error = function(e) {
          warning(paste('Failed to create poly.2:',e))
        })
      }
    }
  }
  
  plt <- ggplot2::ggplot(data=coords, ggplot2::aes(x=X,y=Y, col='original')) + ggplot2::geom_path() + 
    ggplot2::geom_path(data =line.data, ggplot2::aes(x=x, y=y, col='cutoff')) +
    ggplot2::theme_classic()
  
  if(ratio.first.poly > 0 && ratio.first.poly<1) {
    plt <- plt + ggplot2::geom_path(data =line.intersect, ggplot2::aes(x=X, y=Y, col='cutoff after cut'))
  }
  
  if(!is.null(poly.1)){
    plt <- plt + ggplot2::geom_path(data =sf::st_coordinates(poly.1), ggplot2::aes(x=X, y=Y, col='poly.1'))
  }
  
  if(!is.null(poly.2)){
    plt <- plt + ggplot2::geom_path(data =sf::st_coordinates(poly.2), ggplot2::aes(x=X, y=Y, col='poly.2'))
  }
  
  return(list(
    ratio.first.poly = ratio.first.poly,
    full.poly = poly,
    poly.1 = poly.1,
    poly.2 = poly.2,
    plot = plt
  ))
}

#' @title plot_population
#'  
#' @description
#' 
#'    Plot a simulated population. Plot has polygons drawn as lines, with the individuals from population as red dots if is.random=FALSE, blue triangles if is.random=TRUE.
#'   
#' @param polygons list of sf polygons - drawn as lines on the plot
#' @param population results dataframe from any of the generating functions (clumped, random, uniform). Has columns X, Y, and is.random (is the point of a random subpopulation)
#' @param pct_random Percent of the population that is random - for the plot title
#' @param curr_n Number of individuals - for the plot title
#' @param dot_size Default=2, size of the dots on the plot
#' @param wanted.aspect.ratio Default=1, wanted aspect ratio of the plot
#' @param draw.poly Default=TRUE, should polygon be drawn around points
#' @param with.text Default=TRUE, should text (title, legend, axes) appear or not
#' @param coloramp optional colormap as a named vector with random, clumped, and uniform. If NULL, default is used
#'                 c(random='blue', clumped='red', uniform='darkgreen')
#'
#' @return the ggplot2 plot
#' @export
plot_population <- function(polygons, population, pct_random, curr_n, dot_size=2, wanted.aspect.ratio=1,
                            draw.poly=TRUE, with.text = TRUE, colormap = NULL) {
  if(is.null(colormap)) {
    colormap <- c(random='blue', clumped='red', uniform='darkgreen', subrandom='black')
  }
  coords <- sf::st_coordinates(polygons[[1]])
  if(with.text){
    plt <- ggplot2::ggplot(coords, ggplot2::aes(x=X,y=Y)) 
      if(draw.poly){
        plt <- plt + ggplot2::geom_path()
      }
      plt <- plt + ggplot2::geom_point(data=population, ggplot2::aes(x=X,y=Y,colour = population.type), size=dot_size) +
      ggplot2::scale_color_manual(values=colormap) +
      ggplot2::labs(title= paste('Simulated pop - ', pct_random,'% clumped, n=',curr_n)) + 
      ggplot2::theme_classic() + ggplot2::theme(aspect.ratio = wanted.aspect.ratio)
  } else {
      plt <- ggplot2::ggplot(coords, ggplot2::aes(x=X,y=Y))
      if(draw.poly){
        plt <- plt + ggplot2::geom_path()
      }
      plt <- plt + 
        ggplot2::geom_point(data=population, ggplot2::aes(x=X,y=Y,colour = population.type), size=dot_size) +
        ggplot2::scale_color_manual(values=colormap) +
        ggplot2::theme_classic() + ggplot2::theme(aspect.ratio = wanted.aspect.ratio,
                                                  axis.title.x = ggplot2::element_blank(),
                                                  axis.text.x = ggplot2::element_blank(),
                                                  axis.line.x = ggplot2::element_blank(),
                                                  axis.title.y = ggplot2::element_blank(),
                                                  axis.text.y = ggplot2::element_blank(),
                                                  axis.line.y = ggplot2::element_blank(),
                                                  axis.ticks = ggplot2::element_blank(),
                                                  legend.position = "none",
                                                  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "mm"))
  }
  return(plt)
}

#' @title generate_random_population
#'  
#' @description
#' 
#'   Generate a random population of size n within a list of sf polygons.
#' 
#' @param n number of individuals in the population
#' @param polys list of sf polygons in which to randomize the population
#' @return A dataframe of the random population, with columns TAG (numeric, 1-n), X, and Y.
#' 
#' @export
generate_random_population <- function(n, polys) {
  coords <- sf::st_coordinates(polys[[1]])
  xrange <- c(min(coords[,1]),max(coords[,1]))
  yrange <- c(min(coords[,2]),max(coords[,2]))
  individual.locations.data.frame.ret <- data.frame(TAG= 1:n, X=0, Y=0)
  for(j in 1:length(individual.locations.data.frame.ret$X)){
    poly.ind <- sample(1:length(polys),1)
    poly <- polys[[poly.ind]]
    x <- sample((xrange[1]:xrange[2])*10,1)/10
    y <- sample((yrange[1]:yrange[2])*10,1)/10
    while(length(sf::st_contains(poly, sfheaders::sf_point(c(x,y)))[[1]])!=1){
      x <- sample((xrange[1]:xrange[2])*10,1)/10
      y <- sample((yrange[1]:yrange[2])*10,1)/10
    }
    individual.locations.data.frame.ret[j,]$X <- x
    individual.locations.data.frame.ret[j,]$Y <- y
  }
  
  individual.locations.data.frame.ret$population.type <- population_type$random

  return(individual.locations.data.frame.ret)
}


#' @title generate_clumped_population
#'  
#' @description
#' 
#'   Generate a clumped population of n points within polygons in polys, return as a dataframe. Currently the resulting population is split into clumps of sqrt(n) points.
#' 
#' @param n number of individuals in the population
#' @param polygons list of sf polygons in which to randomize the population
#' @param n_for_clumping used for calculating the expected distance for the clumping. Can cause closer clumpedness with n_for_clumping > n
#' @param n_clumps number of clumps wanted in the clumped populations. Default=0 uses the default function of floor(sqrt(n))
#' @param max_dist_from_centroid_divisor divisor of the expected distance for distances of each point from it's randomized centroid. 
#'                                   Higher number causes smaller distances between points in the same clump.
#' @param min_dist_between_clumps_multiplier the multiplier of the same distance for the distances between the centroids/clumps.
#' @return A dataframe of the clumped population, with columns TAG (numeric, 1-n), X, Y, and clump.id (running index of the clump of the individual).
#' 
#' @export
generate_clumped_population <- function(n, polygons, 
                                         n_for_clumping=0, n_clumps=0,
                                         max_dist_from_centroid_divisor=2,
                                         min_dist_between_clumps_multiplier=1.5){
  
  coords <- sf::st_coordinates(polygons[[1]])
  xrange <- c(min(coords[,1]),max(coords[,1]))
  yrange <- c(min(coords[,2]),max(coords[,2]))
  
  w <- xrange[2] - xrange[1]
  h <- yrange[2] - yrange[1]
  
  area.mcp.rect <- w*h
  
  ## for simplicity - use first polygon for now.
  poly <- polygons[[1]]
  
  area.mult <- area.mcp.rect / sf::st_area(poly)
  
  n_for_d <- n
  if(n_for_clumping>0) {
    n_for_d <- n_for_clumping
  }
  
  ##d <- (area.mcp.rect / ceiling(area.mult*n))^ 0.5
  d <- (sf::st_area(poly) / n_for_d)^ 0.5
  
  re <- list()
  individual.locations.data.frame.ret <- data.frame(TAG=1:n, X=0, Y=0)
  individual.locations.data.frame.ret$clump.id <- 0
  
  n_centroids <- n_clumps
  if(n_centroids == 0){
    n_centroids <-  floor(sqrt(n))
  }
  
  j_used <- 1
  
  clump.centroids<-as.data.frame(matrix(nrow=0,ncol=2))
  colnames(clump.centroids) <- c('X','Y')
  
  for(i in 1:n_centroids){
    
    ## RANDOMIZE centroid
    re_randomize <- TRUE
    attempts<-1
    while(re_randomize && attempts<100){
      x_centroid <- sample((xrange[1]:xrange[2])*10,1)/10
      y_centroid <- sample((yrange[1]:yrange[2])*10,1)/10
      
      re_randomize <- length(sf::st_contains(poly, sfheaders::sf_point(c(x_centroid,y_centroid)))[[1]])!=1
      if(i>1 && !re_randomize){
        relevant.centroids <- clump.centroids[1:(i-1),]
        relevant.centroids$dist.from.curr <- ((relevant.centroids$X - x_centroid)^2 + (relevant.centroids$Y - y_centroid)^2)^0.5
        if(length(relevant.centroids[relevant.centroids$dist.from.curr < min_dist_between_clumps_multiplier*d,]$X)>0){
          re_randomize <- TRUE
        }
      }
      attempts<-attempts+1
    }
    
    clump.centroids <- rbind(clump.centroids,c(x_centroid, y_centroid))
    colnames(clump.centroids) <- c('X','Y')
    
    
    for(j in (i-1)*(n%/%n_centroids) + 1:(n%/%n_centroids)){
      
      if(j>n) break
      a <- sample(((1:360)-1)*pi/180,1)
      x <- clump.centroids[i,]$X + sample(1:(d%/%max_dist_from_centroid_divisor),1)* cos(a)
      y <- clump.centroids[i,]$Y + sample(1:(d%/%max_dist_from_centroid_divisor),1)* sin(a)
      attempts <- 1
      while(length(sf::st_contains(poly, sfheaders::sf_point(c(x,y)))[[1]])!=1 && attempts<100){
        a <- sample(((1:360)-1)*pi/180,1)
        x <- clump.centroids[i,]$X + sample(1:(d%/%max_dist_from_centroid_divisor),1)* cos(a)
        y <- clump.centroids[i,]$Y + sample(1:(d%/%max_dist_from_centroid_divisor),1)* sin(a)
        attempts <- attempts +1
        if(attempts==100){
          stop(paste('Too many attempts to generate clumped population with',max_dist_from_centroid_divisor,'- consider dividing by 2'))
          max_dist_from_centroid_divisor <- max_dist_from_centroid_divisor/2
          attempts <-1
        }
      }
      individual.locations.data.frame.ret[j,]$X <- x
      individual.locations.data.frame.ret[j,]$Y <- y
      individual.locations.data.frame.ret[j,]$clump.id <- i
      j_used<-j
    }
  }
  
  while(j_used<n)
  {
    for(i in 1:n_centroids){
      if(j_used==n) break
      j_wanted <- j_used+1
      a <- sample(((1:360)-1)*pi/180,1)
      x <- clump.centroids[i,]$X + sample(1:(d%/%max_dist_from_centroid_divisor),1)* cos(a)
      y <- clump.centroids[i,]$Y + sample(1:(d%/%max_dist_from_centroid_divisor),1)* sin(a)
      while(length(sf::st_contains(poly, sfheaders::sf_point(c(x,y)))[[1]])!=1){
        a <- sample(((1:360)-1)*pi/180,1)
        x <- clump.centroids[i,]$X + sample(1:(d%/%max_dist_from_centroid_divisor),1)* cos(a)
        y <- clump.centroids[i,]$Y + sample(1:(d%/%max_dist_from_centroid_divisor),1)* sin(a)        
      }
      individual.locations.data.frame.ret[j_wanted,]$X <- x
      individual.locations.data.frame.ret[j_wanted,]$Y <- y
      individual.locations.data.frame.ret[j_wanted,]$clump.id <- i
      j_used<-j_wanted
    }
  }
  
  individual.locations.data.frame.ret$population.type <- population_type$clumped
  
  return(individual.locations.data.frame.ret)
}


#' @title generate_uniform_population
#'
#' @description
#' 
#'   Generate a uniform/regular population of n points within polygons in polys, return as a dataframe.
#'      Distance between individuals will be \deqn{UniformDist = \sqrt{\frac{AreaOfPolygons}{n}} \pm noise.ratio * \sqrt{\frac{AreaOfPolygons}{n}}  }
#' 
#' @param n number of individuals in the population
#' @param polygons list of sf polygons in which to randomize the population
#' @param noise.ratio Default 0.1, distance of noise (multiplied by expected distance UniformDist above) to add to each point (so it wont be an exact checkerboard)
#' @param max_dist_from_centroid_divisor divisor of the expected distance for distances of each point from it's randomized centroid. 
#'                                   Higher number causes smaller distances between points in the same clump.
#' @param min_dist_between_clumps_multiplier the multiplier of the same distance for the distances between the centroids/clumps.
#' @return A dataframe of the clumped population, with columns TAG (numeric, 1-n), X, Y, and clump.id (running index of the clump of the individual).
#' 
#' @export
generate_uniform_population <- function(n, polygons, n.for.dist=0, noise.ratio=0.1){
  
  ## for simplicity - use first polygon for now.
  poly <- polygons[[1]]
  
  coords <- sf::st_coordinates(polygons[[1]])
  xrange <- c(min(coords[,1]),max(coords[,1]))
  yrange <- c(min(coords[,2]),max(coords[,2]))
  
  w <- xrange[2] - xrange[1]
  h <- yrange[2] - yrange[1]
  
  ##d <- (area.mcp.rect / ceiling(area.mult*n))^ 0.5
  d_x <- (sf::st_area(poly) / ifelse(n.for.dist==0,n,n.for.dist))^ 0.5
  d_y <- d_x
  noise <- noise.ratio*min(d_x, d_y)
  
  xlim <- max(floor(w/d_x),0)
  ylim <- floor(h/d_y)
  
  if(h<w)
  {
    xlim <- floor(w/d_x)
    ylim <- max(floor(h/d_y),0)
  }
  
  # fix if not correctly divisible
  if(xlim*ylim < n) {
    #print(paste('fix needed - ', xlim, ylim, xlim*ylim, n, w, h))
    missing <- n - xlim*ylim
    
    if(xlim<ylim){
      xlim <- xlim + ceiling(missing/ylim)
    } else {
      ylim <- ylim + ceiling(missing/xlim)
    }
    
    # if(missing > max(ylim, xlim)){
    #     xlim <- xlim + missing%/%ylim
    #     ylim <- ylim + missing%/%xlim
    # } else {
    #  if(xlim<ylim){
    #    xlim <- xlim + missing/ylim
    #  } else {
    #    ylim <- ylim + missing/xlim
    #  }
    # }

    
    d_x <- w/xlim
    d_y <- h/ylim
    noise <- noise.ratio*min(d_x, d_y)
    
    #print(paste('after fix - ', xlim, ylim, xlim*ylim, n, w, h))
    
  }
  
  
  re <- list()
  individual.locations.data.frame.ret <- data.frame(TAG=1:n, X=0, Y=0)
  
  j <- 1
  for(x in xrange[1]-d_x/2+(d_x*(1:xlim))){
    if(j>n) break
    for(y in yrange[1]-d_y/2+(d_y*(1:ylim))){
      if(j>n) break
      x_rand <- x+sample((-noise:noise)*10,1)/10
      y_rand <- y+sample((-noise:noise)*10,1)/10
      attempts<-1
      pip <- length(sf::st_contains(poly, sfheaders::sf_point(c(x_rand,y_rand)))[[1]])==1
      while(attempts<100 && !pip){
        x_rand <- x+sample((-noise:noise)*10,1)/10
        y_rand <- y+sample((-noise:noise)*10,1)/10
        pip <- length(sf::st_contains(poly, sfheaders::sf_point(c(x,y)))[[1]])==1
        attempts <- attempts + 1
      }
      
      if(pip){
        individual.locations.data.frame.ret[j,]$X <- x_rand
        individual.locations.data.frame.ret[j,]$Y <- y_rand
        j <- j+1
      }
    }
  }
  
  if(j<n+1)
  {
    stop(paste('WARNING - Uniform randomization yielded',j-1,'points instead of',n,'\n'))
  }
  
  individual.locations.data.frame.ret$population.type <- population_type$uniform
  
  return(individual.locations.data.frame.ret)
}


#' @title simulate_population_remove_from_original_for_smaller_pops
#'
#' @description
#' 
#' Core function of the process, does  both the heatmap simulation data and the null-model distribution data (via the wrapper function below).
#' The general flow is as follows:
#'        \itemize{
#'        \item{} {Generate simulations_per_percent populations of size initial_n of each level of initial_pct/percent_random_rise + 1, such that each such
#'            simulated population is X% clumped and (100-X)% random. Calculate all the functions in function_names on the populations.}
#'        \item{} {Until the populations dwindle to minimum size 2, remove delt_n individuals from each population at each step, and recalculate
#'            all the functions in function_names.}
#'        }
#' This is meant to test the robustness of the functions in function_names of differentiating between different levels of
#' clumpedness for even small subpopulations, where it is assumed that only some of the real population have been geotagged.
#' Note that this function parallely runs on N-1 out of N cores of the PC, and that the heavy runtime parameter is set by a mix of 
#' percent_random_rise and simulations_per_percent, since these set the original number of populations that need to be simulated.
#' 
#' @param polygons list of sf polygons where the populations can be generated.
#' @param function_names vector of function names (Clumping Methods, CMs), which are first tested by check_function_validity.
#' @param helper_functions_to_export vector of additional functions needed for the CMs to run. Tested only for existence in memory.
#' @param initial_n number of individual of each original population Initial Sample Size (ISS) in the paper.
#' @param n.for.calc sent into each of the function in function_names, but does not have to be used by each of these. Is meant
##                to normalize for the tested initial_n, since the null-model distribution data is suggested to be used on population
##                sizes which may be larger than initial_n. For example - if a quadrat.test is used, the null-model distribution
##                calculation MUST take into account the amount of quadrats used for the smaller initial_n.
#' @param n.for.clumping similarly, the clumping function above calculates distances between centroids and between individuals
##                according to the number of individuals. So in this same logic, using the higher number of individuals in the 
##                null-model distribution will cause much smaller distances within clumps. To normalize if wanted, this can be set
##                to a non-zero value - a zero value causes the number of individuals sent to the clumping function to be used. The 
##                change of using and not using this seems minor in any case.
#' @param n.for.uniform default = 0, if 0 then the CURRENT amount is taken (ie - the uniform subpopulation is spread over the entire polygon),
#'                      if nonzero the value itself is taken.
#' @param sim.type default = simulation_type$random.to.clumped, which part of the range to simulate.
#' @param uniform.random.generation.type default = uniform_random_generation_type$split.polygon.and.generate.separately, how to generate the populations in the uniform-random range (DFR -100 to 0), if relevant.
#' @param uniform.noise.ratio default = 0.1, noise for uniform function
#' @param initial_div divisor of expected distance of individuals from the centroid of their clump in the clumped sub-populations.
#' @param n_clumps number of clumps wanted in the clumped populations. Default=0 uses the default function of floor(sqrt(n))
#' @param delt_n the amount of individuals to remove from each of the simulated population for each step. Higher amount lowers sensitivity,
##              so it will be harder to ascertain at which number of individuals exactly the measures start to fail. Removed Individuals (RIs) in the paper.
#' @param initial_pct the maximum percentage of clumping to start from - for direct usage of this function this would normally be 100%,
##              but for example the null-model distribution would choose low clumpedness percentages to check, ie a value of 20 maximum.
#' @param percent_random_rise the drops in percent that set the levels of clumpedness - like above, this sets the sensitivity of the
##              second axis, ie how good the measures are at differentiating between different clumpedness levels. However - this also
##              incurrs a higher price on runtime, since it also sets the amount of initial populations that need to be generated. Sets the levels of Degree of Clumping (dfr) from the paper.
#' @param simulations_per_percent the amount of simulated populations generated at each level of clumpedness at initial_n, and therefore
##              also for each removal of delt_n individuals. Simulation Repetitions (SRs) in the paper.
#' @param recalculate_poly default=T, should the polygons be recalculated to be the polygon of the remaining population. 
#'              Use this to encapsulate tagging with an unknown real potential population.
#' @param recalculate_poly_with_rect default=F, should the recalcualted polygon (if T for above) should be a bounding rectangle, or the convex polygon.
#' @param plotResults as the name implies, whether to plot and display one simulated population for each level of clumpedness using ggplot.
#' @param ns_to_plot which ns should be plotted. If empty, will be only the original population size (ISS)
#' @param n_cores Default=0, number of cores to use for parallel processing. A value of zero will use the benchmark of parallel:detectCores(logical=FALSE)-1.
#'                Results of benchmark tests with initial_n = 30, simulations_per_percent = 100, delt_n = 5, and percent_random_rise = 20 
#'                are below - all values over 7 cores took less than a minute.
#'                
#'                  2 cores \tab 1.158628 mins \tab 1.431186 mins \tab 1.447352 mins\cr
#'                  3 cores \tab 1.450939 mins \tab 1.615264 mins \tab 1.551022 mins\cr
#'                  4 cores \tab 1.267707 mins \tab 1.110228 mins \tab 1.097605 mins\cr
#'                  5 cores \tab 1.016974 mins \tab 1.0203 mins \tab 1.021301 mins\cr
#'                  6 cores \tab 1.29498 mins \tab 1.368751 mins \tab 1.404914 mins\cr
#'                  7 cores \tab 1.263218 mins \tab 1.3216 mins \tab 59.86292 secs\cr
#'                  8 cores \tab 59.03853 secs \tab 59.46413 secs \tab 57.9628 secs\cr
#'                  9 cores \tab 59.18239 secs \tab 57.93962 secs \tab 58.88627 secs\cr
#'                  10 cores \tab 57.72106 secs \tab 1.039918 mins \tab 57.46983 secs\cr
#'                  11 cores \tab 57.00588 secs \tab 57.55166 secs \tab 58.5118 secs\cr
#'                  12 cores \tab 57.31099 secs \tab 1.042125 mins \tab 55.91966 secs\cr
#'                  13 cores \tab 56.57559 secs \tab 57.63825 secs \tab 1.015002 mins\cr
#'                  14 cores \tab 58.69229 secs \tab 57.11787 secs \tab 58.8794 secs\cr
#'                  15 cores \tab 59.30693 secs \tab 58.96996 secs \tab 59.86175 secs\cr
#'                  16 cores \tab 59.50359 secs \tab 59.90103 secs \tab 59.55054 secs
#' 
#' @return a dataframe with length(function_names) + 2 columns and 
#'  \deqn{numRows = simulationsPerPercent \times {\frac{initialN}{deltN}} \times (1+{\frac{initialPct}{percentRandomRise}})}
#'         rows, such that the first two columns css and dfr are the number of individuals in the population and the percent clumpedness of one population,
#'         and the remaining columns hold the output of the function function_names_i for the population.
#' 
#' @export
simulate_population_remove_from_original_for_smaller_pops <- function(polygons, 
                                                                      function_names,
                                                                      helper_functions_to_export,
                                                                      sim.type = simulation_type$random.to.clumped,
                                                                      uniform.random.generation.type = uniform_random_generation_type$split.polygon.and.generate.separately,
                                                                      initial_n = 20,
                                                                      n.for.calc = 0,
                                                                      n.for.clumping = 0,
                                                                      n.for.uniform = 0,
                                                                      initial_div = 2,
                                                                      n_clumps = 0,
                                                                      uniform.noise.ratio = 0.1,
                                                                      delt_n=5,
                                                                      initial_pct=100,
                                                                      percent_random_rise=20,
                                                                      simulations_per_percent=20,
                                                                      recalculate_poly=T,
                                                                      recalculate_poly_with_rect=F,
                                                                      plotResults=FALSE,
                                                                      ns_to_plot = c(),
                                                                      n_cores = 0) {
  function_list <- c()
  
  failures <- character()
  
  for(x in function_names){
    ret.check <- check_clumping_method_validity(polygons, x, sim.type)
    if(!ret.check$success){
      failures <- paste0(failures,'\n * ',ret.check$message)
    } else {
      function_list <- c(function_list,get(x))
    }
  }
  
  for(y in helper_functions_to_export) {
    func <- NULL
    
    ## check if function exists
    tryCatch({
      func <- get(y)
    }, error = function(e) {
      failures <- paste0(failures,'\n * ', paste('Helper Function',y,'does not exist in environment'))
    })
    
    ## check if null at this point sometimes get does not collapse
    if(is.null(func)){
      failures <- paste0(failures,'\n * ', paste('Helper Function',y,'does not exist in environment'))
    }
  }
  
  if(length(failures)>0){
    failures <- paste0('The following failures happened:', failures)
    stop(failures)
  }
  
  final = matrix(nrow=0,ncol=2 + length(function_names))
  
  curr_n = initial_n
  curr_pct_random = 0
  
  ns_to_plot = ns_to_plot
  if(length(ns_to_plot)==0){
    ns_to_plot <- c(initial_n)
  }
  
  ## local copy for R - for parallel.
  polygons <- polygons
  n.for.calc <- n.for.calc
  n.for.clumping <- n.for.clumping
  
  if(n_cores==0){
    n_cores <- parallel::detectCores(logical = F)-1
  }
  
  cluster <- parallel::makeCluster(n_cores)
  
  parallel::clusterExport(cluster,c(c('generate_simulation_polygons','generate_clumped_population', 
                                      'generate_random_population', 'generate_uniform_population',
                                      'population_type','simulation_type', 'uniform_random_generation_type'),
                            helper_functions_to_export,
                            function_names))
  
  initial_pops <- NULL
  
  dfr.vector <- c()
  
  if(sim.type!=simulation_type$uniform.to.random.to.clumped) {
    dfr.vector <- 0:(simulations_per_percent*(1+initial_pct/percent_random_rise)-1)
  } else {
    dfr.vector <- 0:(simulations_per_percent*(1+2*initial_pct/percent_random_rise)-1)
  }
  
  split.polys <- list()
  
  need.split.polys <- FALSE
  
  if(uniform.random.generation.type == uniform_random_generation_type$split.polygon.and.generate.separately &&
    sim.type != simulation_type$random.to.clumped) 
  {
      need.split.polys <- TRUE
      num.splits.needed <- (1+initial_pct/percent_random_rise)
      print(paste('Note - will generate',num.splits.needed,'split polygons, since sim.type=',sim.type,'and uniform.random.generation.type requested demands split.'))
      
      for(i in 1:num.splits.needed) {
        curr.ratio.uniform <- ((i-1)*percent_random_rise)/100
        split.polys <- append(split.polys,list(split_poly(polys[[1]],curr.ratio.uniform)))
      }
  }
  
  ## min possible for CHULL is 4
  while(curr_n >= delt_n && curr_n > 4){
    
    if(!is.null(initial_pops)){
      gs <- parallel::parLapply(cluster,dfr.vector, function(i) {
        curr_pct <- percent_random_rise*(i%/%simulations_per_percent)
        curr_id <- i%%simulations_per_percent
        
        if(sim.type==simulation_type$uniform.to.random){
          curr_pct <- curr_pct * -1
        } else if (sim.type==simulation_type$uniform.to.random.to.clumped) {
          curr_pct <- curr_pct - initial_pct
        }
        
        inds <- sample(1:(curr_n+delt_n),curr_n)
        
        ret <- list()
        ret$id <- curr_id
        ret$pct_random <- curr_pct
        ret$curr.population <- initial_pops[[i+1]]$curr.population[inds,]
        
        if(recalculate_poly) {
          polygons.for.funcs <- generate_simulation_polygons(ret$curr.population, use_bounding_rect = recalculate_poly_with_rect)
        } else {
          polygons.for.funcs <- polygons
        }
        
        ret$results <- c()
        for(f in function_list){
          res <- f(ret$curr.population, polygons.for.funcs, n.for.calc)
          ret$results <- c(ret$results, res)
        }
        return(ret)
      })
      
      initial_pops <- gs
    }
    else {
      gs <- parallel::parLapply(cluster,dfr.vector, function(i) {
        ret <- list()
        clumped.subpopulation <- NULL
        random.subpopulation <- NULL
        
        curr_pct <- percent_random_rise*(i%/%simulations_per_percent)
        curr_id <- i%%simulations_per_percent
        
        if(sim.type==simulation_type$uniform.to.random){
          curr_pct <- curr_pct * -1
        } else if (sim.type==simulation_type$uniform.to.random.to.clumped) {
          curr_pct <- curr_pct - initial_pct
        }

        amt_nonrandom <- floor(curr_n * abs(curr_pct)/100)
        
        ## clumped
        
        nonrandom.subpopulation <- NULL
        
        if(need.split.polys && curr_pct<0) {
          split.poly <- split.polys[lapply(split.polys,function(x) {x$ratio.first.poly==abs(curr_pct)/100})==TRUE][[1]]
          if(amt_nonrandom>0){
            nonrandom.subpopulation <- generate_uniform_population(amt_nonrandom, list(split.poly$poly.1), noise.ratio = uniform.noise.ratio)
            nonrandom.subpopulation$clump.id <- -1
          }
          
          ## random
          if(amt_nonrandom<initial_n){
            random.subpopulation <- generate_random_population(initial_n - amt_nonrandom, list(split.poly$poly.2))
            random.subpopulation$clump.id <- 0
          }
        } else {
            if(amt_nonrandom>0){
              
              if(curr_pct>=0) {
                    nonrandom.subpopulation <- generate_clumped_population(amt_nonrandom, polygons, n_clumps = n_clumps,
                                                      n_for_clumping = ifelse(n.for.clumping==0, initial_n, n.for.clumping), 
                                                      max_dist_from_centroid_divisor = initial_div)
              } else {
                  nonrandom.subpopulation <- generate_uniform_population(amt_nonrandom, polygons, uniform.noise.ratio,
                                                                         n.for.dist = ifelse(n.for.uniform==0, amt_nonrandom, n.for.uniform))
                  nonrandom.subpopulation$clump.id <- -1
              }
            }
            
            ## random
            if(amt_nonrandom<initial_n){
              random.subpopulation <- generate_random_population(initial_n - amt_nonrandom, polygons)
              random.subpopulation$clump.id <- 0
            }
        }
        
        curr.population <- nonrandom.subpopulation
        
        if(is.null(curr.population)){
          curr.population <- random.subpopulation
        } else if (!is.null(random.subpopulation)){
          curr.population <- rbind(curr.population, random.subpopulation)
        }
        
        ret$id <- curr_id
        ret$curr.population <- curr.population
        ret$pct_random <- curr_pct
        ret$results <- c()
        
        if(recalculate_poly) {
          polygons.for.funcs <- generate_simulation_polygons(ret$curr.population, use_bounding_rect = recalculate_poly_with_rect)
        } else {
          polygons.for.funcs <- polygons
        }
        
        for(f in function_list){
          res <- f(curr.population, polygons.for.funcs, n.for.calc)
          ret$results <- c(ret$results, res)
        }
        
        return(ret)
      })
      
      initial_pops <- gs
    }
    
    for(x in gs){
      final <- rbind(final,c(c(curr_n, x$pct_random),x$results))
    }
    
    if(plotResults && curr_n %in% ns_to_plot){
      for(i in 1:length(initial_pops)){
        if(initial_pops[[i]]$id==0){
          plot_population(polygons, initial_pops[[i]]$curr.population, initial_pops[[i]]$pct_random, curr_n)
        }
      }
    }
    
    curr_n <- curr_n - delt_n
  }
  
  parallel::stopCluster(cluster)
  
  
  final.heatmap.data <- as.data.frame(final)
  
  colnames(final.heatmap.data) <- c(c('css','dfr'), function_names)
  
  return(final.heatmap.data)
  
}

#' @title generate_large_randomlike_population_distribution
#'
#' @description
#' 
#' Using simulate_population_remove_from_original_for_smaller_pops, generate a Large Random-Like Population Distribution (LRPD),
#' which is later compared to values from the same function. The reasoning here is that:
#' \itemize{
#' \item{}{A large enough population (default n=100) such that any basic Clumping Method (CM) will be able to differentiate between a very clumped and very random population}
#' \item{}{A low Degree of Clumping (dfr) - such that these large simulated populations should come out on the random end (hence "Random-Like"). Default is 0-20% dfr.}
#' }
#' 
#' 
#' @param polygons list of sf polygons where the populations can be generated.
#' @param function_names vector of function names (Clumping Methods, CMs), which are first tested by check_function_validity.
#' @param helper_functions_to_export vector of additional functions needed for the CMs to run. Tested only for existence in memory.
#' @param n.for.calc sent into each of the function in function_names, but does not have to be used by each of these. Is meant
##                to normalize for the tested initial_n, since the null-model distribution data is suggested to be used on population
##                sizes which may be larger than initial_n. For example - if a quadrat.test is used, the null-model distribution
##                calculation MUST take into account the amount of quadrats used for the smaller initial_n.
#' @param n.for.clumping similarly, the clumping function above calculates distances between centroids and between individuals
##                according to the number of individuals. So in this same logic, using the higher number of individuals in the 
##                null-model distribution will cause much smaller distances within clumps. To normalize if wanted, this can be set
##                to a non-zero value - a zero value causes the number of individuals sent to the clumping function to be used. The 
##                change of using and not using this seems minor in any case.
#' @param n.for.uniform default = 0, if 0 then the CURRENT amount is taken (ie - the uniform subpopulation is spread over the entire polygon),
#'                      if nonzero the value itself is taken.
#' @param sim.type default = simulation_type$random.to.clumped, which part of the range to simulate.
#' @param uniform.random.generation.type default = uniform_random_generation_type$split.polygon.and.generate.separately, how to generate the populations in the uniform-random range (DFR -100 to 0), if relevant.
#' @param uniform.noise.ratio default = 0.1, noise for uniform function
#' @param initial_div divisor of expected distance of individuals from the centroid of their clump in the clumped sub-populations.
#' @param delt_n the amount of individuals to remove from each of the simulated population for each step. Higher amount lowers sensitivity,
##              so it will be harder to ascertain at which number of individuals exactly the measures start to fail. Removed Individuals (RIs) in the paper.
#' @param max_pct_to_use the maximum percentage of clumping to start from - for direct usage of this function this would normally be 100%,
##              but for example the null-model distribution would choose low clumpedness percentages to check, ie a value of 20 maximum.
#' @param percent_random_rise the drops in percent that set the levels of clumpedness - like above, this sets the sensitivity of the
##              second axis, ie how good the measures are at differentiating between different clumpedness levels. However - this also
##              incurrs a higher price on runtime, since it also sets the amount of initial populations that need to be generated. Sets the levels of Degree of Clumping (dfr) from the paper.
#' @param simulations_per_percent the amount of simulated populations generated at each level of clumpedness at initial_n, and therefore
##              also for each removal of delt_n individuals. Simulation Repetitions (SRs) in the paper.
#' @param recalculate_poly default=T, should the polygons be recalculated to be the polygon of the remaining population. 
#'              Use this to encapsulate tagging with an unknown real potential population.
#' @param recalculate_poly_with_rect default=F, should the recalcualted polygon (if T for above) should be a bounding rectangle, or the convex polygon.
#' @param plotResults as the name implies, whether to plot and display one simulated population for each level of clumpedness using ggplot.
#' @param n_cores Default=0, number of cores to use for parallel processing. A value of zero will use the benchmark of parallel:detectCores(logical=FALSE)-1.
#' 
#' @return a dataframe with length(function_names) + 2 columns and 
#'  \deqn{numRows = simulationsPerPercent \times (1+{\frac{initialPct}{percentRandomRise}})}
#'         rows, such that the first two columns css and dfr are the number of individuals in the population and the percent clumpedness of one population,
#'         and the remaining columns hold the output of the function function_names_i for the population.
#' 
#' @export
generate_large_randomlike_population_distribution <- function(polygons,
                                                              function_names,
                                                              helper_functions_to_export,
                                                              sim.type = simulation_type$random.to.clumped,
                                                              uniform.random.generation.type = uniform_random_generation_type$split.polygon.and.generate.separately,
                                                              uniform.noise.ratio = 0.1,
                                                              n.for.calc,
                                                              n.for.clumping=0,
                                                              n.for.uniform=0,
                                                              initial_div=2,
                                                              max_pct_to_use=20, 
                                                              percent_random_rise = 10,
                                                              simulations_per_percent=100,
                                                              recalculate_poly=T,
                                                              recalculate_poly_with_rect=F,
                                                              plotResults = F,
                                                              n_cores = 6) {
  
  df.res.distrib <- simulate_population_remove_from_original_for_smaller_pops(polygons, 
                                                                              function_names,
                                                                              helper_functions_to_export,
                                                                              sim.type = sim.type,
                                                                              uniform.random.generation.type = uniform.random.generation.type,
                                                                              uniform.noise.ratio = uniform.noise.ratio,
                                                                              n.for.calc=n.for.calc,
                                                                              n.for.clumping=n.for.clumping,
                                                                              n.for.uniform=n.for.uniform,
                                                                              initial_n = 100, delt_n = 99, ## 100 and then stop
                                                                              plotResults = plotResults, 
                                                                              simulations_per_percent = simulations_per_percent,
                                                                              initial_div = initial_div, 
                                                                              percent_random_rise = percent_random_rise,
                                                                              initial_pct = max_pct_to_use,
                                                                              recalculate_poly = recalculate_poly,
                                                                              recalculate_poly_with_rect = recalculate_poly_with_rect,
                                                                              n_cores = n_cores)
  return(df.res.distrib)
}

#' @title calculate_nonrandom_percent_on_heatmap_data
#'
#' @description
#' For each Clumping Method (CM) of columns 3 and above in both heatmap.data and lrpd.data:
#' \itemize{
#' \item{}{Calculate the Minimum Non-Random Threshold (MNRT) - the upperPercentile percentile of the lrpd.data for the CM}
#' \item{}{Append a boolean column to the heatmap.data, such that 1 means a CM value over the MNRT}
#' }
#' The returned value will contain the Non-Random Percent (NRP) per Current Sample Size (CSS) and Degree of Clumping (dfr) pair.
#' 
#' @param heatmap.data A result dataframe of simulate_population_remove_from_original_for_smaller_pops
#' @param lrpd.data A result dataframe of generate_large_randomlike_population_distribution - containing a Large Random-Like Population Distribution for each CM
#' @param upperPercentile Default=0.95 (95th percentile), the percentile of the LRPD to be taken as the MNRT per CM. For negative DRF, 1-upperPercentile will be taken as maximum instead.
#' @param d_cutoff Which dfr (Degree of Clumping) percentages to keep for the final MNRT
#' 
#' @return An aggregated (mean) dataframe, with the columns of heatmap.data with and additional column for every CM in heatmap.data (column 3 and above), such that the name of the column
#'          for a CM named 'name' will be called 'name.non.random.percent'. The results are the mean value (in the original column) and the mean
#'          Non-Random Percent (in the new columns) per CSS + dfr pair (Current Sample Size and Degree of Clumping respectively).
#'          This data can later be plotted in various ways to test the robustness of each CM.
#' 
#' @export
calculate_nonrandom_percent_on_heatmap_data <- function(heatmap.data,
                                                      lrpd.data,
                                                      upperPercentile=0.95,
                                                      d_cutoff = c(0))
{
  for(col in colnames(heatmap.data)){
    if(!col %in% c('css','dfr')){
        heatmap.data[paste0(col,'.non.random.percent')] <- -1
        num.pos <- nrow(lrpd.data[lrpd.data$dfr>0,])
        num.neg <- nrow(lrpd.data[lrpd.data$dfr<0,])
        
        if(num.pos>0){
            mnrt <- quantile(lrpd.data[lrpd.data$dfr %in% d_cutoff & lrpd.data$dfr>=0,col],c(upperPercentile))
            heatmap.data[heatmap.data$dfr>0,paste0(col,'.non.random.percent')] <- as.numeric(heatmap.data[heatmap.data$dfr>0,col]>mnrt)
            if(num.neg==0){
              heatmap.data[heatmap.data$dfr==0,paste0(col,'.non.random.percent')] <- as.numeric(heatmap.data[heatmap.data$dfr==0,col]>mnrt)
            }
        }
        if(num.neg>0){
            d_cutoff <- d_cutoff * -1
            mnrt <- quantile(lrpd.data[lrpd.data$dfr %in% d_cutoff & lrpd.data$dfr<=0,c(col)],c(1-upperPercentile))
            heatmap.data[heatmap.data$dfr<0,paste0(col,'.non.random.percent')] <- as.numeric(heatmap.data[heatmap.data$dfr<0,col]<mnrt)
            heatmap.data[heatmap.data$dfr==0,paste0(col,'.non.random.percent')] <- as.numeric(heatmap.data[heatmap.data$dfr==0,col]<mnrt)
        }
    }
  }
  
  agg <- aggregate(heatmap.data,by=list(heatmap.data$css, heatmap.data$dfr),FUN=function(x) {mean(x)})
  agg <- agg[,3:ncol(agg)]
  
  return(agg)
}

#' @title calculate_values_example_data
#'
#' @description
#' On example data bird.location.data, for each population, calculate the grade of each Clumping Method (CM) given by funcs.
#' 
#' @param example_individual_location_data intended for example data bird.location.data, but can be used for any dataframe with the following columns: POPULATION (character), X, Y.
#' @param funcs vector of function names. Will be tested for validity as a Clumping Method.
#' @param check_clumping_method_valid default=TRUE, check if clumping method is valid or not (via check_clumping_method_validity)
#' 
#' @return A dataframe with the following columns: POPULATION, css (amount of individuals), and a column with the value of every CM on the population.
#'
#' @export
calculate_values_example_data <- function(example_individual_location_data, funcs, check_clumping_method_valid=T){
  
  pops <- unique(example_individual_location_data$POPULATION)
  return_data <- data.frame(POPULATION=pops,css=0, area=0)
  for(f in funcs){
    return_data[f] <- NaN
  }
  
  for(i in 1:length(pops)){
    sub <- example_individual_location_data[example_individual_location_data$POPULATION == pops[i], ]
    n <- nrow(sub)
    return_data[i,c('css')] <- n
    polys <- generate_simulation_polygons(sub,use_bounding_rect = F)
    return_data[i,c('area')] <- sf::st_area(polys[[1]]) 
    
    if(i==1 && check_clumping_method_valid){
      failures <- character()
      
      for(x in funcs){
        ret.check <- check_clumping_method_validity(polys, x)
        if(!ret.check$success){
          failures <- paste0(failures,'\n * ',ret.check$message)
        }
      }
      
      if(length(failures)>0){
        failures <- paste0('The following failures happened:', failures)
        stop(failures)
      }
    }
    
    for(f in funcs){
      fun <- get(f)
      return_data[i,c(f)] <- fun(sub, polys, n)
    }
  }
  return(return_data)
}
