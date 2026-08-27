


library(ggplot2)

library(tidyr)

#source("scripts/Uniformity_tests.R")

#library(bayesplot)

library(latex2exp)

library(glue)









# close form shapley values for an average value function v.



shapley_mean_closedform <- function(x) {

  n <- length(x)

  if (n == 0) return(numeric(0))

  Hn <- sum(1 / seq_len(n))

  Sh <- numeric(n)

  for (i in seq_len(n)) {

    mean_others <- if (n==1) 0 else sum(x[-i]) / (n - 1)

    Sh[i] <- (1 / n) * x[i] + ((Hn - 1) / n) * (x[i] - mean_others)

  }

  return(Sh)

}





cauchy_space_piet <- function(x){  

  pe <- pexp(-log(x), rate = 1)

  ps <- 2*pmin(pe,1-pe)

 return( tan((0.5-ps)*pi) )

}



cauchy_space_pot <- function(x){

  pb <- pbeta(sort(x), 1:length(x), seq(length(x), 1, by=-1))

  ps <-2*pmin(pb, 1-pb)

  return(  tan((0.5-ps)*pi) )



  }





cauchy_space_prit<- function(x) {

  N <- length(x)             

  ranks <-colSums(outer(x, x, "<=")) # R=N.F(x) : evaluate and scale ecdf at each point x_i 

  probs1<-pbinom(ranks,N,x)          # P(R <= r_i) ,  R follows Binom(N,x)

  probs2 <-pbinom(ranks-1,N,x)       # P(R>r_i-1)=1-P(R <= r_i-1)

  p_vals <-2*pmin(probs1, 1-probs2)  # individual test pvalues        

  tan((0.5-p_vals)*pi)                

 }







# Cauchy combination

cauchy_agg <- function(x) { 1 - pcauchy(mean(tan((0.5-x)*pi))) } # x vector ( dependent pits)



# Truncate Cauchy combination

Tcauchy_agg <- function(x) {

  mask <- as.numeric(x < 0.5)  # truncation mask

  1 - pcauchy(mean(tan((0.5-x)*pi)*mask) ) }





influential_points_idx<- function(x,alpha=0.05){



target <- qcauchy(1-alpha)



# sort positive values descending and keep their original indices

pos_idx <- order(x, decreasing = TRUE)

pos_idx <- pos_idx[x[pos_idx] > 0]

pos_vals <- x[pos_idx]



# cumulative removal sums

cumsum_remove <- cumsum(pos_vals)



# find minimal number of removals to reach target

needed <- which(sum(x) - cumsum_remove <= target)[1]



# indices and values to remove

removed_idx <- pos_idx[seq_len(needed)]

return(removed_idx)

#removed_vals <- x[removed_idx]

#  remaining vector

  #remaining <- x[-removed_idx]

}



##########################################################################################







################  POT-C & PIET-C & PRIT-C  tests ######################





# POT-C : Pointwise Order-based tests Combination ( Beta-based tests)

pot_c<- function(x) {

 probs <- pbeta(sort(x), 1:length(x), seq(length(x), 1, by=-1))

 ps <-2*pmin(probs, 1-probs)

 Tcauchy_agg(ps) }





# PRIT-C : Pointwise Rank-based Individual Tests Combination ( Binomial-based tests)

prit_c<- function(x) {

   N <- length(x)             

  ranks <-colSums(outer(x, x, "<=")) # R=N.F(x) : evaluate and scale ecdf at each points x_i 

  probs1<-pbinom(ranks,N,x)          # P(R <= r_i) ,  R follow Binom(N,x)

  probs2 <-pbinom(ranks-1,N,x)       # P(R>r_i-1)=1-P(R <= r_i-1)

  p_vals <-2*pmin(probs1, 1-probs2)  # individual test p-values        

  Tcauchy_agg(p_vals)                # p-values aggregation

 }







## PIET-C : Pointwise Inverse-CDF Evaluation Tests Combination ( Exp(1)-based tests)

piet_c<- function(x) {

  pe <- pexp(-log(x), rate = 1)

  ps <- 2*pmin(pe,1-pe)

  cauchy_agg(ps)}




## ECDF-based Visual representation of POT-C,  PRIT-C and PIET-C.



 graphical_rep<-function(x,test="POT-C", xlab="LOO-PIT", gamma=NULL,

                         infl_points_only=FALSE,diff=FALSE,alpha=0.05){



  if (test == "POT-C") {

    cauchy_vals <- cauchy_space_pot(sort(x))

    pval <- pot_c(x)

  } else if (test == "PIET-C") {

    cauchy_vals <- cauchy_space_piet(sort(x))

    pval  <-piet_c(x)

  } else if (test=="PRIT-C"){

    cauchy_vals <- cauchy_space_prit(sort(x))

    pval  <-prit_c(x)

  } else {stop("Unknown test: available options are test='POT-C', test='PRIT-C' and test='PIET-C'" )

  }



   # calculate expected marginal contributions of u_(i) to the overall test.

  sh_val <- shapley_mean_closedform(cauchy_vals)



  

   color <- ifelse(pval<0.05, "red", "#F97316")

    

   if(diff==TRUE){

     ECDF <-colMeans(outer(x, x, "<="))-x

   }else{ECDF <- colMeans(outer(x, x, "<="))}

   

  df <- data.frame(

    pit = x,

    ecdf =ECDF # colMeans(outer(x, x, "<="))-x

  )



  # Sort by pit

  df <- df[order(df$pit), ]



  p <- ggplot(df, aes(x = pit, y = ecdf)) +

    geom_line(color = "black", linewidth = 0.5) +          

    geom_point(color = "black", size = 0.1) +          



 # geom_abline(intercept = 0, slope = 0,

             # linetype = "solid", color = "darkgray", linewidth = 0.7 ,alpha=0.6) +

    labs(               #test,

      title =paste(paste("Uniformity p-value ="),floor(pval * 1000) / 1000 ) ,

      x = xlab,

      y = "ECDF"

    ) +

  theme_minimal(base_size = 16)+

    theme(

     panel.grid = element_blank(),

    axis.line = element_line(color = "black"), # show black axis lines

    axis.ticks = element_line(color = "black"),# show ticks

    axis.text = element_text(color = "black",size=16), # keep axis labels visible

    axis.title = element_text(color = "black",size=16),

    plot.title = element_text(size = 15,face = "bold"),

    )


   if(diff==TRUE){ p <- p+ geom_abline(intercept = 0, slope = 0,

                linetype = "solid", color = "darkgray", linewidth = 0.7 ,alpha=0.6)+

                scale_y_continuous(limits=c(-0.1,0.1))+

                labs(y = "ECDF-difference")+

                scale_x_continuous(limits = c(0,1.03),breaks = c(0,0.25,0.50,0.75,1),expand = c(0.01,0))

                  }else{p <- p+geom_abline(intercept = 0, slope = 1,

                           linetype = "solid", color = "darkgray", linewidth = 0.7 ,alpha=0.6) +

                           scale_x_continuous(limits=c(0,1.03),breaks = c(0,0.25,0.50,0.75,1))+

                           scale_y_continuous(limits=c(0,1.03),breaks = c(0.25,0.50,0.75,1))+

                           coord_cartesian(expand = FALSE) }

                                       

  if(infl_points_only==TRUE){

    if(pval<0.05){

      pos_idx <- influential_points_idx(x=sh_val,alpha=alpha)

      df_points <- df[pos_idx, ]

      p <- p+ geom_point( data = df_points, color = color,size = 1.4)

    message(

  sprintf(

    "%d highly influential values detected in the ordered values at positions : %s",

    length(pos_idx),

    paste(pos_idx, collapse = ", ")

  )

)



    }

   #if(pval>=0.05){p <- p}

  } else

     {

 ## Identify suspecious points (or regions)

 

        #tolerence

  if (is.null(gamma) || gamma<0){

  gamma <- ifelse(pval <= 0.5, 0,

                  ifelse(pval <= 0.9, 0.05, 0.2)) }

  

  red_idx <- which(sh_val>gamma)

       

  if(length(red_idx)!=0){ 

    

    df_red <- df[red_idx, ]



    # groups of consecutive  suspecious point 

    consec_groups <- cumsum(c(1, diff(red_idx) != 1))

    df_red$segment <- consec_groups   # assign group label

    

    # Separate isolated vs grouped points

    df_isolated   <- df_red[ ave(df_red$pit, df_red$segment, FUN = length) == 1, ]

    df_grouped    <- df_red[ ave(df_red$pit, df_red$segment, FUN = length) > 1, ]
  

    p <-p + geom_line( data = df_grouped, aes(group = segment),color = color,linewidth =1.5) +

      # Red dots ONLY on isolated points

      geom_point( data = df_isolated, color = color,size = 1.7)



    

    message(sprintf("tolerance level (gamma), \u03b3 \u2208 [0, %s] : \u03b3=%s",

                    floor(max(sh_val) * 100) / 100,gamma )

            )

}



  }

     

return(p) # return(pval) 

}



##################################################################################################





############# Examples 







## exemple 1 : evidence against uniform ( rejected) : problematic regions in red

## set.seed(2049)



## x <- 1 - (1 - runif(300))^(1.2)



## ## POT-C 

## pot_c(x) # test p-value



## x11()

## graphical_rep(x,test="POT-C", infl_points_only=FALSE)












