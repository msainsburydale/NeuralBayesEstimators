suppressMessages({
library("ggplot2")
library("ggpubr")
library("dplyr")
})

combine_xi <- function(xi_train, xi_test) {
  train_df <- xi_train %>% data.frame %>% cbind(Parameters = "Training")
  test_df  <- xi_test %>% data.frame %>% cbind(Parameters = "Testing")
  
  return(rbind(train_df, test_df))
}

plot_theta_lambda <-  function(xi_train, xi_test) {
  df <- combine_xi(xi_train, xi_test)
  ggplot(df, aes(x = logLambda, y = exp(logTheta))) + 
    geom_point(aes(colour = Parameters), size = 0.1, alpha = 0.5) + 
    labs(x = expression(log(lambda)), y = expression(theta)) + 
    theme_bw()
}


plot_theta_nu <-  function(xi_train, xi_test) {
  df <- combine_xi(xi_train, xi_test)
  ggplot(df, aes(x = logNu, y = logTheta)) + 
    geom_point(aes(colour = Parameters), size = 0.7, alpha = 0.5) + 
    labs(x = expression(log(nu)), y = expression(log(theta))) + 
    theme_bw()
}

plot_theta_lambda_nu <- function(y, separate_panels = FALSE, use_all_params = FALSE) {
  
  y <- as.data.frame(y)
  
  plot_of_lambda_and_nu <- y %>% 
    ggplot() + 
    geom_point(aes(x = logTheta, y = logNu)) + 
    labs(x = expression(log(theta)), y = expression(log(nu))#, 
         # title = as.expression(bquote(
         #   "All combinations of " ~ theta ~ "and" ~ log(nu)
         # ))
         ) + 
    theme_bw() #+ 
    #theme(plot.title = element_text(hjust = 0.5)) 
  
  
  ## Split by the smoothness parameter
  theta_lambda_by_nu <- split(y, y$logNu)
  
  if (use_all_params) {
    idx <- 1:length(theta_lambda_by_nu)
  } else {
    ## Just plot the first, (rough) middle, and last value of nu:
    idx <- floor(seq(1, length(theta_lambda_by_nu), length.out = 3))
  }

  
  y_subset <- theta_lambda_by_nu[idx] %>% 
    abind(along = 1) %>% 
    as.data.frame
  
  if (separate_panels) {
    ## Compute the range of lambda so that we can enforce the same x-axis, which 
    ## will facilitate comparison between the panels:
    lambda_range <- y_subset %>% 
      select("logLambda") %>% 
      range
    
    plotlist <- lapply(theta_lambda_by_nu[idx], function(x) {
      ggplot(x, aes(x = logLambda, y = logTheta)) +
        geom_point(size = 0.1, alpha = 0.5) +
        labs(x = expression(log(lambda)), y = expression(log(theta)),
             title = as.expression(bquote(
               nu ~ "=" ~ .(round(exp(x$logNu[1]), 3))
             ))) +
        xlim(lambda_range) +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5)) 
    })
  } else {
    
    if (!use_all_params) {
      ## Convert to a factor for prettier plotting
      y_subset$logNu <- y_subset$logNu %>% 
        round(3) %>% 
        as.factor
    }

    plotlist <- ggplot(y_subset, aes(x = logLambda, y = logTheta, colour = logNu)) +
      geom_point(size = 0.1, alpha = 0.5) +
      labs(x = expression(log(lambda)), y = expression(log(theta)), 
           colour = expression(log(nu))) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5)) 
    
    if (use_all_params) plotlist <- plotlist + scale_colour_distiller(palette = "Spectral")
    
    plotlist <- list(plotlist)
  }
  
  
  
  plotlist <- c(list(plot_of_lambda_and_nu), plotlist)
  
  figure <- ggpubr::ggarrange(plotlist = plotlist, nrow = 1, align = "hv")
  
  return(figure)
}


## Original code taken from: 
## https://github.com/bertcarnell/lhs/blob/master/vignettes/lhs_basics.Rmd

graph2dLHS <- function(Alhs)
{
  stopifnot(ncol(Alhs) == 2)
  sims <- nrow(Alhs)
  par(mar = c(4,4,2,2))
  plot.default(Alhs[,1], Alhs[,2], type = "n", ylim = c(0,1),
               xlim = c(0,1), xlab = "Parameter 1", ylab = "Parameter 2", xaxs = "i", 
               yaxs = "i", main = "")
  for (i in 1:nrow(Alhs))
  {
    rect(floor(Alhs[i,1]*sims)/sims, floor(Alhs[i,2]*sims)/sims,
         ceiling(Alhs[i,1]*sims)/sims, ceiling(Alhs[i,2]*sims)/sims, col = "grey")
  }
  points(Alhs[,1], Alhs[,2], pch = 20, col = "red")
  abline(v = (0:sims)/sims, h = (0:sims)/sims)
}

# transform is a function of the kind that takes a number
# transform <- function(x){return(qnorm(x,mean=0, std=1))}
graph2dLHSTransform <- function(Alhs, transform1, transform2, 
                                min1 = NULL, max1 = NULL, min2 = NULL, max2 = NULL, 
                                drawgrid = TRUE)
{
  if (is.null(min1)) min1 <- min(Alhs[, 1])
  if (is.null(max1)) max1 <- max(Alhs[, 1])
  if (is.null(min2)) min2 <- min(Alhs[, 2])
  if (is.null(max2)) max2 <- max(Alhs[, 2])
  
  stopifnot(ncol(Alhs) == 2)
  stopifnot(all(Alhs[,1] <= max1 && Alhs[,1] >= min1))
  stopifnot(all(Alhs[,2] <= max2 && Alhs[,2] >= min2))
  sims <- nrow(Alhs)
  breaks <- seq(0,1,length = sims + 1)[2:(sims)]
  breaksTransformed1 <- sapply(breaks, transform1)
  breaksTransformed2 <- sapply(breaks, transform2)
  
  par(mar = c(4,4,2,2))
  plot.default(Alhs[,1], Alhs[,2], type = "n", 
               ylim = c(min2, max2),
               xlim = c(min1, max1),
               xlab = "Parameter 1", ylab = "Parameter 2", 
               xaxs = "i", yaxs = "i", main = "")
  
  if (drawgrid) {
    for (si in 1:sims)
    {
      temp <- Alhs[si,]
      for (i in 1:sims)
      {
        if ((i == 1 && min1 <= temp[1] && breaksTransformed1[i] >= temp[1]) ||
            (i == sims && max1 >= temp[1] && breaksTransformed1[i - 1] <= temp[1]) ||
            (breaksTransformed1[i - 1] <= temp[1] && breaksTransformed1[i] >= temp[1]))
        {
          for (j in 1:sims)
          {
            if ((j == 1 && min2 <= temp[2] && breaksTransformed2[j] >= temp[2]) ||
                (j == sims && max2 >= temp[2] && breaksTransformed2[j - 1] <= temp[2]) ||
                (breaksTransformed2[j - 1] <= temp[2] && breaksTransformed2[j] >= temp[2]))
            {
              if (i == 1)
              {
                xbot <- min1
                xtop <- breaksTransformed1[i]
              } else if (i == sims)
              {
                xbot <- breaksTransformed1[i - 1]
                xtop <- max1
              } else 
              {
                xbot <- breaksTransformed1[i - 1]
                xtop <- breaksTransformed1[i]
              }
              if (j == 1)
              {
                ybot <- min2
                ytop <- breaksTransformed2[j]
              } else if (j == sims)
              {
                ybot <- breaksTransformed2[j - 1]
                ytop <- max2
              } else 
              {
                ybot <- breaksTransformed2[j - 1]
                ytop <- breaksTransformed2[j]
              }
              rect(xbot, ybot, xtop, ytop, col = "grey")
            }
            
          }
        }
      }
    }
    
    abline(v = breaksTransformed1, h = breaksTransformed2)
  }
  
  col <- if (!drawgrid) scales::alpha("red", 0.2) else "red"
  
  points(Alhs[,1], Alhs[,2], pch = 20, col = col)
}

#set.seed(1111)
#A <- randomLHS(5,4)
#f <- function(x){qnorm(x)}
#g <- function(x){qlnorm(x, meanlog=0.5, sdlog=1)}
#B <- A
#B[,1] <- f(A[,1])
#B[,2] <- g(A[,2])
#graph2dLHSTransform(B[,1:2], f, g, -4, 4, 0, 8)
#f <- function(x){qunif(x, 3, 5)}
#B <- apply(A, 2, f)
#graph2dLHSTransform(B[,1:2], f)
