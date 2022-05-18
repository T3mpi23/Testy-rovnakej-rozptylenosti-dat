
# Moodova testová štatistika
Mood_statistic <- function(X,Y)
{
  
  W <- 0
  m <- length(X)
  n <- length(Y)
  
  for (i in (m+1):(m+n))
  {
    W <- W + (  (  rank( c(X,Y) )[i]  ) - (m+n+1)/2 )^2
  }
  
  return (W)
}

# Ansariho-Bradleyho testová štatistika
Ansari_statistic <- function(X,Y)
{
  rank_vector <-c()
  even <- 0
  m <- length(X)
  n <- length(Y)
  s <- m+n
  C <-0
  
   rank_vector <- rank(c(X,Y))  #spojená sad dát
  
  if (s %%2==0) 
  {
     
    #pridelenie rankov
    rank_vector[rank_vector > s/2 ] <- (s+1)- rank_vector[rank_vector > s/2 ]
     
    E_C <- n*(s+2)/4
    D_C <- m*n*(s+2)*(s-2)/(48*(s-1))
    
    }
  if (s %%2==1) 
  {
   
    #pridelenie rankov
    rank_vector[rank_vector>(s+1)/2 ]<-(s+1)-rank_vector[rank_vector >(s+1)/2]
    
    E_C <- n*(s+1)^2 / (4*s)
    D_C <- m*n*(s+1)*(s^2 +2)/(48*(s^2))
  }

  C<- sum( rank_vector[(m+1):s])
  return(C)
  
}

#Lepageova testová štatistika
Lepage_statistic <- function(X,Y )
{
  m <- length(X)
  n <- length(Y)
  
  #Wilcoxonová testová štatistika
  W = sum(rank(c(Y,X))[1:n] ) 
  E_W <-n*(m+n+1)/2  #stedná hodnota E(W)
  D_W <- m*n*(m+n+1)/12 #rozptyl D(W)

  s <- m+n
   
  #Stredná hodnota a rozptyl Ansariho-Bradleyho test. štat.
  if (s %%2==0) 
  {
    E_C <- n*(s+2)/4
    D_C <- m*n*(s+2)*(s-2)/(48*(s-1))
    
  }
  if (s %%2==1) 
  {
    E_C <- n*(s+1)^2 / (4*s)
    D_C <- m*n*(s+1)*(s^2 +2)/(48*(s^2))
  }
  C <- Ansari_statistic(X,Y)
  #Lepageova test. štatistika
  D <- ((W -E_W)^2 )/D_W  + ((C-E_C)^2)/D_C
  return(D)
}

#Millerova testová štatistika
Miller_statistic <- function(X,Y )
{
  m <- length(X)
  n <- length(Y)
  
  #pomocné matice, kde každý riadok je sadá dát X (resp. Y)
  #Odstránením diagonály získame sady dát v ktorých postupne 
  #vynechávame prvky (Jacknife)
  X_matrix <- matrix( rep(X,m), byrow = TRUE, ncol=m)
  Y_matrix <- matrix( rep(Y,n), byrow = TRUE, ncol=n)
  
     diag(X_matrix)<-NA
  X_matrix<-t(matrix(t( X_matrix)[which(!is.na( X_matrix))],nrow=m-1,ncol=m))
  
  diag(Y_matrix)<-NA
  Y_matrix<-t(matrix(t( Y_matrix)[which(!is.na( Y_matrix))],nrow=n-1,ncol=n))
 
   library("matrixStats")
  V_X <- rowVars(X_matrix)
  V_Y <- rowVars(Y_matrix)
   E_X <- rowMeans(X_matrix)
   E_Y <- rowMeans(Y_matrix)
  S_i <- log(V_X)
  T_i <- log(V_Y)
  S_0 <- log(var(X))
  T_0 <- log( var(Y))

  A <-  m* rep(S_0, m) - (m-1)*S_i 
  B <- n*rep(T_0,n)  - (n - 1)*T_i
  A_mean <- mean(A)
  B_mean <- mean(B)
    V_1 <- var(A)/m
  V_2 <- var(B)/n
  
  Q <- (A_mean - B_mean)/ sqrt(V_1 + V_2)
  
return(Q)   
}
#Lehmanova testová štatistika 
Lehman_statistic <-function(X,Y){
  
   m<- length(X)
   n<- length(Y)
   
   # Kroneckerov produkt pre urýchlenie, urobí všetky rozdiely
    X_diffs <-  matrix ( abs(kronecker(X,X,FUN ="-")), 
                         byrow = TRUE, nrow = m)
    X_diffs<- X_diffs[ upper.tri(X_diffs)]
    
    Y_diffs  <-  matrix ( (abs(kronecker(Y,Y,FUN ="-")  )), 
                          byrow = TRUE, nrow = n)
    Y_diffs <- Y_diffs[upper.tri(Y_diffs)]
    # rozdiely rozdielov, berieme len pozitivné
   M <- kronecker(X_diffs,Y_diffs, FUN = "-" )[kronecker(X_diffs,
                                                         Y_diffs, FUN = "-")>0] 
   #Lehmannova test. štat.
   L <- length(   M   )
   
   
   
   return(L)
}
# Permutaèný prístup odhadu disperzie Lehmannovej test. štat.
# B je poèet kôl
# výstup: vektor Lehmanových test. štat.
# odhad disperzie je potom rozptyl výstupu
Permutaions <- function(X,Y,B)
{
  L_perm <-c()
  X_rand <-c()
  Y_rand <-c()
  m <- length(X)
  n <- length(Y)
  
  for (b in  1:B)
  {
    rand_index <- sample(1:( m+n), length(X)+length(Y), replace =  FALSE)
    X_rand <- c(X,Y)[rand_index[1:m]]
    Y_rand <- c(X,Y)[rand_index[(m+1):(n+m)]]
    L_perm[length(L_perm)+1] <- Lehman_statistic(X_rand, Y_rand)
  }
  return(L_perm)
}
# Bootstrap odhad Lehmannovej test. štat.
#B poèet kôl
# výstup: vektor Lehmanových test. štat.
bootstrap  <- function(X,Y,B) 
  {
  L_boot <-c()
  X_rand <-c()
  Y_rand <-c()
  m <- length(X)
  n <- length(Y)
  for (b in  1:B)
  {
    # tvorba nových sád dát
    rand_index <- sample(1:( m+n), length(X)+length(Y), replace =  TRUE)
      X_rand <- c(X,Y)[rand_index[1:m]]
      Y_rand <- c(X,Y)[rand_index[(m+1):(n+m)]]
    L_boot[length(L_boot)+1] <- Lehman_statistic(X_rand, Y_rand)
  }
  return( (L_boot) )
}
# modifikovaný bootstrap kde medzi sebou nemiešame prvky z rôznych sád
bootstrap_individual <- function(X,Y,B)
{
  L_boot <-c()
  X_rand <-c()
  Y_rand <-c()
  m <- length(X)
  n <- length(Y)
  for (b in  1:B)
  {
    rand_X_index <- sample(1:m, m, replace =  TRUE)
    rand_Y_index <- sample(1:n, n, replace =  TRUE)
    X_rand <- X[rand_X_index[1:m]]
    Y_rand <- Y[rand_Y_index[1:n]]
    L_boot[length(L_boot)+1] <- Lehman_statistic(X_rand, Y_rand)
  }
  
  return( L_boot )
}
# bootstrap ja zlepšenie odhadu strednej hodnoty
bootmean <- function(X,B)
{
  X_rand <-c()
  m <- length(X)
  Mean_boot<-c()
  for (b in  1:B)
  {
    rand_index <- sample(1:m, m, replace =  TRUE)
    X_rand <- X[rand_index[1:m]]
    Mean_boot[length(Mean_boot)+1] <- mean(X_rand)
  }
  return(Mean_boot)
}



#Paralel setup
##########
#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}
parallel::detectCores()
n.cores <- parallel::detectCores() - 2

#loading example data
data("penguins")

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


parallel::stopCluster(cl = my.cluster)

closeparalels <-function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
closeparalels()

#########


#1Testing Type I error
######################
# 
# critical <- qnorm(0.975,mean=0,sd=1)
# chi_critical <- qchisq(0.95, 2)
# #means_norm_dist
#   mx <- 0 
#   my <- 0
# #location of cauchy
#   loc_X <- 0
#   loc_Y <- 0
#   
# #pomocne var.  
# count_F_norm <-0
# count_F_cauchy <-0
# count_M_norm <- 0
# count_M_cauchy <- 0
# 
# m = 100 #size of X data
# n = 100 # size of Y data   Pre ANSARI TREBA ABY n <= m 
# 
# for (i in 1:10000)
# {
#   
#   #a F-test 
#   #Normal dist.
#   X_norm <- rnorm(m,mean = mx, sd = 1)
#   Y_norm <- rnorm(n, mean = my, sd = 1)
#   
#   if (var.test(X_norm,Y_norm, conf.level = 0.95)$p.value < 0.05)
#   {
#     count_F_norm <- count_F_norm + 1
#   }
#   #cauchy dist.
#   X_cauchy <- rcauchy(m, loc_X, 1)
#   Y_cauchy <- rcauchy(n, loc_Y, 1)
#   
#   if (var.test(X_cauchy,Y_cauchy, conf.level = 0.95)$p.value < 0.05)
#   {
#     count_F_cauchy <- count_F_cauchy + 1
#   }
#   
#   
#   #b Mood-test 
#   ##########################xxx
#   s <- m + n
#   #for norm dist.
#   M_norm <- 0
#  
#   for (i in m:(m+n))
#   {
#     M_norm <- M_norm + (  (  rank( c(X_norm,Y_norm) )[i]  ) - (m+n+1)/2 )^2
#   }
#   E_M_norm <- n*(s+1)*(s-1)/12
#   D_M_norm <- m*n*(s+1)*(s+2)*(s-2)/180
#   
#   if ( ( (M_norm-E_M_norm)/sqrt(D_M_norm) > critical ) | ((M_norm-E_M_norm)/sqrt(D_M_norm) < (-1*critical)) ){
#     count_M_norm <- count_M_norm + 1
#   }
#   
#   #for cauchy dist.
#   M_cauchy <- 0
#   
#   for (i in m:(m+n))
#   {
#     M_cauchy <- M_cauchy + (  (  rank( c(X_cauchy,Y_cauchy) )[i]  ) - (m+n+1)/2 )^2
#   }
#   E_M_cauchy <- n*(s+1)*(s-1)/12
#   D_M_cauchy <- m*n*(s+1)*(s+2)*(s-2)/180
#   
#   if ( ( (M_cauchy-E_M_cauchy)/sqrt(D_M_cauchy) > critical ) | ((M_cauchy-E_M_cauchy)/sqrt(D_M_cauchy) < (-1*critical)) ){
#     count_M_cauchy <- count_M_cauchy + 1
#   }
#   
# }
# 
# Type_I_err_F_norm <- count_F_norm/10000
# Type_I_err_F_cauchy <- count_F_cauchy/10000
# Type_I_err_M_norm <- count_M_norm/10000
# Type_I_err_M_cauchy <- count_M_cauchy/10000
# 
# 





#####################################################################################################################################
# powers of methods

  

    

      replic <-1000

          ay_Fisher <- c()
          ay_Mood <-c()
  
          ay_Ansari_asym <-c()
          ay_Ansari_table <-c()
          ay_Ansari_MYtable <- c()
          
          
           ay_Lehman_asym <-c()
  
          ay_Miller <-c()
          ay_Miller_k2 <-c()
          ay_Miller_rand <-c()
          ay_Miller_stud <-c()
  
          ay_Lepage_asym <-c()
          ay_Lepage_table <-c()
  
          ay_Lehman_pivot <-c()
          ay_Lehman_pivot_indiv <-c()
  
          ay_Lehman_perm <-c()
          ay_Lehman_perm_centered <-c()
          ay_Lehman_perm_centered_boot <-c()
  
  
          ay_Lehman_individ_bootstrap <- c()

        m = 30#size of X data
          n = 30# size of Y data   Ansari n <= m   (m + n)<= 20
         
       
           s <- m + n
         
           c_AB_upper <- 68  #rozdelime alfa1 + alfa2 = alfa/2 ,  P(C > c_upper) = 0.0255  reject if more
          c_AB_lower <- 42 - 1
          
           # original_lower <- 26 - 1
           # original_upper <- 51
           
        
           d_lepage_critical <- 5.718615



           student_criticals <- qt(0.975,df=m+n-2)
        
           lehman_criticals_upper <-1471
             lehman_criticals_lower <- 557- 1
             
           critical <- qnorm(0.975,mean=0,sd=1)
        chi_critical <- qchisq(0.95, 2)
        
          mx <-0
        my <-0
        
        sdx_max <-2# ide o disperzie aj ked su oznacene ako sd 
        sdy<-1
        tick<-0.1
        osX<-seq(0.40, sdx_max,by=tick)
        
        
        

      
          boot_rep <-200 #number   of bootstrap reps fo  r pivotal bootsrap
          Student_mat<- matrix(nrow = 5, ncol = 54)
           TIE <- c()
          
           d1 <- Sys.time()
            set.seed(2)   
         
        
             for (sdx in osX)
              {
              count_Fisher <- 0
              count_Mood <- 0
              
              count_Ansari_asym <- 0
              count_Ansari_table <- 0
              count_Ansari_MYtable <- 0
              
              count_Lehman_asym <- 0
              count_Lehman_pivor <- 0 
              count_Miller <- 0
              count_Miller_k2 <-0
              count_Miller_rand <-0
              
              count_Lepage_asym <- 0
              count_Lepage_table <- 0
                 
             # for( i in 1:replic)
              
       
     
                u <- foreach(i = 1:replic, .combine = "cbind") %dopar%
                  {
                #Generating ///////
                # Normal
  # 

                # X <- rnorm(m,mean = mx, sd = sqrt(sdx))
                # Y <- rnorm(n, mean = my, sd = sqrt(1))
              # # cauchy
              
                
                bool_Fisher <-0
                bool_Mood <-0
                bool_Ansari <-0
                bool_Miller <-0
                bool_Lepage <- 0
                bool_Lehman <-0
                  
                
                # X <- rcauchy(m, location = mx,scale = sdx )
                # Y <- rcauchy(n, location = my,scale = 1)
                X <- rt(m, df= 5)*(sqrt(sdx))
                Y <- rt(n, df = 5)
                # X_mean <- X - mean(X)
                # Y_mean <- Y - mean(Y)
                # 
                # X <- X_mean
                # Y <- Y_mean
               
                # X_median <- X - median(X)
                # Y_median <- Y - median(Y)
                # 
                # X <- X_median
                # Y <- Y_median
                
                #FISHER TEST
                #############################
                # if (var.test(X,Y, conf.level = 0.95)$p.value <= 0.05)
                # {
                #   bool_Fisher <-1
                # # count_Fisher <- count_Fisher + 1
                # }
                # ##############################
              
                #Moood block
                #############################

                W <- Mood_statistic(X,Y)
                #
                E_W <-n*(s+1)*(s-1)/12
                D_W <- m*n*(s+1)*(s+2)*(s-2)/180
                if ( ( abs(W-E_W)/sqrt(D_W) >= critical ) )
                  {
                  # count_Mood <- count_Mood + 1
                  bool_Mood <-1
                }
                # 
                #                
                #                   
                                  
                # #end of mood blok
                ##############################
                
                # ANSARI-Bradley blok
                #########################################
                C <- Ansari_statistic(X,Y)
                # # # C_med <- Ansari_statistic(X_median,Y_median)
                if (s %%2==0)
                {
                  E_C <- n*(s+2)/4
                  D_C <- m*n*(s+2)*(s-2)/(48*(s-1))
                }
                if (s %%2==1)
                {
                  E_C <- n*(s+1)^2 / (4*s)
                  D_C <- m*n*(s+1)*(s^2 +3)/(48*(s^2))
                }
#
# #
                  if ( abs( ( C-E_C )/sqrt(D_C) ) >= critical )
                  {
                    bool_Ansari <-1
#                     # count_Ansari_asym <- count_Ansari_asym + 1
                  }

                # if ( (C >= original_upper )  |  ( C <= original_lower) )
                # {
                #   bool_Ansari <-1
                #   # count_Ansari_table <- count_Ansari_table + 1
                # }
# #

                  # if (   (C >= c_AB_upper )  |  ( C <= c_AB_lower)  ){
                  #   bool_Ansari<-1
                  #   # count_Ansari_MYtable <-count_Ansari_MYtable + 1
                  # }


                
                # if (   (C_med >= c_AB_upper )  |  ( C_med <= c_AB_lower)  ){
                #   bool_Fisher <-1
                #   # count_Ansari_MYtable <-count_Ansari_MYtable + 1
                # }
                
                ##########################################
                # LEHMAN
                ######################################################
             
                # L<-Lehman_statistic(X,Y)
                # E_L <- m*(m-1)*n*(n-1)/8
                # D_L_perm <- var (Permutaion_test(X,Y,300))
                # 
                # if(  (  abs ((L-E_L)/sqrt(D_L_perm))  >= critical )  )
                # {
                #   bool_Lehman <- 1
                # }
  
               

                # D_L_perm <- var (Permutaion_test(X,Y,k*200))
                # if(  (  abs ((L-E_L)/sqrt(D_L_perm))  >= critical )  )
                # {
                #   bool_Lehman <- 1
                #   # count_Lehman_perm_centered = count_Lehman_perm_centered + 1
                # }
  
                # 
                # X_center_boot <- X - mean( bootmean(X,1000))
                # Y_center_boot <- Y - mean( bootmean(Y,1000))
                # D_L_perm_center_boot <- var (Permutaion_test(X_center_boot,Y_center_boot,50))
                # if(  (  abs ((L-E_L)/sqrt(D_L_perm_center_boot))  > critical )  )
                # {
                #   count_Lehman_perm_centered_boot = count_Lehman_perm_centered_boot + 1
                # }
                # 
               
                # D_L_boot <- var(bootstrap(X ,Y ,200))
                # if(  (  abs ((L-E_L)/sqrt(D_L_boot))  >= critical )  )
                # {
                #   bool_Lehman <- 1
                # # count_Lehman_asym = count_Lehman_asym + 1
                # }

                
             
                # L_med <- Lehman_statistic(X_median,Y_median)
                # D_L_boot <- var(bootstrap(X ,Y,200))
                # if(  (  abs ((L-E_L)/sqrt(D_L_boot))  >= critical )  )
                # {
                #   bool_Lehman <- 1
                #   # count_Lehman_asym = count_Lehman_asym + 1
                # }
                # # 
                # #  
            
                
                    # D_L_boot_indi <- var ( bootstrap_individual(X,Y,boot_rep))
               
                # if(  (  abs ((L-E_L)/sqrt(D_L_boot_indi))  > critical )  )
                # {
                #   count_Lehman_individ_bootstrap = count_Lehman_individ_bootstrap + 1
                # }
                # 
                # Pivotal Lehman
                # L_boot <- bootstrap(X,Y,500*k)
                # 
                # # #
                # if( L<= quantile(L_boot, probs = 0.025) | L>= quantile(L_boot, probs = 0.975))
                # {
                # # count_Lehman_pivot <- count_Lehman_pivot + 1
                # bool_Lehman <-1
                # }

                #
               #  if (L <= lehman_criticals_lower | L >= lehman_criticals_upper)
               # {
               #    bool_Lehman <- 1
               #  }
                ######################################################
                
                #Miller test
                ###################
                # # 
                Q <- Miller_statistic(X,Y)
                if ( abs(Q)>=critical )
                  {
                  bool_Miller <- 1
                  # count_Miller <- count_Miller + 1
                  # miller_student_bool <-1
                }

             
  
                # 
                # 
                ###################
                  # ### LEPAGE
                # ###################
                D <- Lepage_statistic(Y,X)
                # # #
                # #
                  if (D >= chi_critical){
                    bool_Lepage <- 1

                  }
                
# 
#                   if (D >= d_lepage_critical){
#                     bool_Lepage <- 1
#                     # count_Lepage_table <- count_Lepage_table + 1
#                   }
                #   # # # #
                # 
                # ###################
  
                
               c(bool_Fisher, bool_Mood,  bool_Ansari,  bool_Lepage,bool_Miller,bool_Lehman)
                
               } #end of for
               
             
          
              # ay_Fisher[ length(ay_Fisher)+1]  <-  sum(u[1,])/replic
              ay_Mood[length(ay_Mood)+1] <-sum(u[2,])/replic
              # 
              ay_Ansari_asym[length(ay_Ansari_asym)+1] <-sum(u[3,])/replic
              ay_Lepage_asym[length(ay_Lepage_asym)+1]<- sum(u[4,])/replic
               ay_Miller[length(ay_Miller)+1]<- sum(u[5,])/replic

               
              ay_Lehman_asym[length(ay_Lehman_asym)+1] <- sum(u[6,])/replic
              # ay_Lehman_perm[length(ay_Lehman_perm)+1] <- sum(u[6,])/replic
              # ay_Lehman_perm_centered[length(ay_Lehman_perm_centered)+1] <- count_Lehman_perm_centered/replic
              # ay_Lehman_perm_centered_boot[length(ay_Lehman_perm_centered_boot)+1] <- count_Lehman_perm_centered_boot/replic
              # ay_Lehman_pivot[length(ay_Lehman_pivot)+1] <-sum(u)/replic
              # ay_Miller_k2[length(ay_Miller_k2)+1]<-  count_Miller_k2/replic
              # ay_Miller_stud[ length(ay_Miller_stud)+1] <- sum(u[2,])/replic
              # ay_Miller_rand[length(ay_Miller_rand)+1]<- count_Miller_rand/replic

               # ay_Lepage_table[length(ay_Lepage_table)+1]<- sum(u[4,])/replic
              # ay_Lehman_individ_bootstrap[ length(ay_Lehman_individ_bootstrap)+1]<- count_Lehman_individ_bootstrap/replic
             u <- c()
               }  
           
            d2 <- Sys.time()
            d2 - d1
        
            
             
        
################################### 
 #GRAPHICAL 

  #plot of whole graph   F test
 
           # sigma
 plot(osX, osX+10, col='blue', ,ylim = c(0,0.3), ylab   = 'Sila',lty=2, xlab = expression(paste(frac(sigma[X]^2,sigma[Y]^2)))         )
        
              # lambda
plot(osX, osX+10, col='blue', ,ylim = c(0,0.7), ylab   = 'Sila',lty=2, xlab = expression(paste(frac(lambda[X],lambda[Y])) ) 
       )
           abline(0.05, 0)
              axis(2, at= 0.05,pos =0.35)   
            
            
            
            points(osX, ay_Fisher, col = 'red')
            lines(osX, ay_Fisher, col = 'red')
            
            points(osX, ay_Lehman_asym, col = 'black')
            lines(osX, ay_Lehman_asym, col = 'black')
            
            points(osX, ay_Miller, col = 'blue')
            lines(osX, ay_Miller, col = 'blue')
            
            points(osX, ay_Ansari_asym, col = 'green')
            lines(osX, ay_Ansari_asym, col = 'green')
            
            points(osX, ay_Mood, col = 'orange')
            lines(osX, ay_Mood, col = 'orange')
            
            points(osX, ay_Lepage_asym, col = 'purple')
            lines(osX, ay_Lepage_asym, col = 'purple')
            
            
         
            
          lines(osX,ay_Ansari_MYtable, col = 'red')
          lines(osX, tmp, col = "red")
          
          lines(osX,ay_Miller, col = 'orange')
          points(osX,ay_Miller, col = 'orange')
          
          lines(osX,ay_Mood, col = 'orange')
          lines(osX,ay_Lehman_asym, col = 'red') 
          lines(osX,ay_Fisher, col = 'green')
          legend(x=0.6,y=0.6,
                 legend=c( 
                          
                          "Lehmannov test",
                          "Moodov test",
                          "Ansariho-Bradleyho test",
                          "Lepageov test"),
                 lty =c(1), col=c("black","orange","green","purple"))
                 
            
          

legend(x=0.6,y=0.3,
       legend=c( "200 kôl bootstrapu","600 kôl bootstrapu","1000 kôl bootstrapu"),
       lty=c(1,1,1,1,1,1),
       col=c("orange","red","blue") )
 
 
plot( (1:22)*500, Lehma_pivot_TIE_cauchy[1:22],ylim = c(0.05,0.1),col="red", xlab = "poèet kôl pri pivotálnom bootstrape"
      , ylab = "odhadnutá pravdepodobnos chyby prvého druhu")
   lines((1:22)*500, Lehma_pivot_TIE_cauchy[1:22], col="red")
 
   n =15
  mx=0
   my = 0
  
  sdx = 2
  boot_rep <- 200
   
 count<-0
  for (i in 1:5000)
  {
    
   X <- rnorm(m,mean = mx, sd = sqrt(sdx))
   Y <- rnorm(n, mean = my, sd = 1)
     A <- sort (bootstrap(X,Y,boot_rep))
     
      a<- Lehman_statistic(X,Y)
      if( a< A[boot_rep*0.025] | a>A[boot_rep *0.975 ] ) count <- count + 1
  }
  
  

   
   
   # Ansari criticals
  #########################
   m = 12
   n = 8
   Lep<-c()
   for (m in 2:12)
   {
     Lep<-c()
     nam <- paste("A", m, sep = "")
     for (n in m:(25-m))
     {
       a<-foreach (i = 1:100000, .combine = 'c') %dopar%
           {
             X<-rnorm(m,0,1)
             Y<-rnorm(n,0,1)
             C <- Lepage_statistic(X,Y)
             C
           }
       Lep <- cbind(Lep,a)
     }
     assign(nam, Lep)
   
   }
 
   
   dim(Ansari_stats)
   quantile(Ansari_stats[,2],probs = 0.025)
   header <- c(6:15)
   Left_axis <- c(3:12)
   
  
    Ansari_critical_0025 <- matrix(ncol = length(header),nrow= length(Left_axis))
    row.names(Ansari_critical_0025) <-Left_axis
    colnames(Ansari_critical_0025) <- header
    Ansari_critical_0025 <- cbind("n", Ansari_critical_0025)
    Ansari_critical_0025 <- rbind("m", Ansari_critical_0025)
    a<- c()
    shift = 0
    ite <-1
   for (n in 2:10)
   {
     for ( n in 3:min(m,25-m))
     {
         Ansari_critical_0025[n-3+1+1,m-4] <-      
        ite <-ite + 1 
     }
    shift <- shift + 1     
   }
   
    
    Ansari_critical_0975 <- matrix(ncol = length(header),nrow= length(Left_axis))
    row.names(Ansari_critical_0975) <-Left_axis
    colnames(Ansari_critical_0975) <- header
    Ansari_critical_0975 <- cbind("n", Ansari_critical_0975)
    Ansari_critical_0975 <- rbind("m", Ansari_critical_0975)
    a<- c()
    shift = 0
    ite <-1
    
    for (m in 6:15)
    {
      for ( n in 3:min(m,25-m))
      {
        Ansari_critical_0975[n-3+1+1,m-4] <- quantile(Ansari_stats[,ite+shift], probs = 0.975)     
        ite <- ite + 1
      }
      shift <- shift + 1     
    }
    
    Ansari_critical_0025_cauchy <- Ansari_critical_0025
      Ansari_critical_0975_cauchy <- Ansari_critical_0975 

      
      
   
      get( paste("A",u,sep = "") )      
     #make matrices, empty
      for (index in 3:12)
      {
        l <- dim(get( paste("A",index,sep = "") )  )[2]
        A <-get( paste("A",index,sep = "") )
        nam <-paste("A_crit_n",index,sep = "")
        tmp <- matrix(ncol = l, nrow = length( min(A):max(A) ) ) 
        # row.names( tmp ) <- paste("x=", c(min(A):max(A)))
        assign(nam,  tmp ) 
      }
      
      #fill them
     for (index in 3:12)
     {
       l <- dim(get( paste("A",index,sep = "") )  )[2]
       A <-get( paste("A",index,sep = "") )
       tmp <- matrix(ncol = l, nrow = length( min(A):max(A) ) ) 
       nam <-paste("Ansari_upper_tail_prob_n",index,sep = "")
       row.names( tmp ) <- paste("x=", c(min(A):max(A)))
       # colnames(tmp) <-paste("m=",c(3:min(index,25-index)))
       colnames(tmp) <-paste("m=",c(index:(25-index)))
       for (c in 1:l)
       {
         vect <- A[,c]
         iter <- 1
         for (crit_val in min(vect):(   max( vect) )   )
         {
               per <-length( vect[vect>=crit_val])/length(vect)
               tmp[iter, c] <- per
               iter<- iter + 1
           
         }
       }
       assign(nam,  tmp )
     }
  
    
      colnames(A3)<- paste("m=",c(3:(25-3)))
      colnames(A4)<- paste("m=",c(4:(25-4)))
      colnames(A5)<- paste("m=",c(5:(25-5)))
      colnames(A6)<- paste("m=",c(6:(25-6)))
      colnames(A7)<- paste("m=",c(7:(25-7)))
      colnames(A8)<- paste("m=",c(8:(25-8)))
      colnames(A9)<- paste("m=",c(9:(25-9)))
      colnames(A10)<- paste("m=",c(10:(25-10)))
      colnames(A11)<- paste("m=",c(11:(25-11)))
      colnames(A12)<- paste("m=",c(12:(25-12)))
      
     
      
      
      
count <- 0
      
  
  for ( n in 3:12)
  {
    nam <-paste("Ansari_upper_tail_prob_n",n,sep = "")
     A <- get( nam)
     A[A <0.01] <- NA
     B <-A[rowSums(!(is.na(A)))!=0,]
     
    assign(nam, B)
  }
  
save( file="Ansari_upper_probablity_tables.RData", Ansari_upper_tail_prob_n3,
            Ansari_upper_tail_prob_n4,
            Ansari_upper_tail_prob_n5, Ansari_upper_tail_prob_n6, Ansari_upper_tail_prob_n7,
            Ansari_upper_tail_prob_n8, Ansari_upper_tail_prob_n9, Ansari_upper_tail_prob_n10,
            Ansari_upper_tail_prob_n11,Ansari_upper_tail_prob_n11,Ansari_upper_tail_prob_n12)



Ansari_upper_tail_prob_final <- matrix()
for ( n in 3:12)
{
  nam <-paste("Ansari_upper_tail_prob_n",n,sep = "")
  A <- get( nam)
  write.csv2( A, file= paste(nam,'.csv',sep = ''))
}


      
  #########################
  
  
  ####### Lepage raw data for critics table 
  ####################
    Lep<-c()
  set.seed(1)
  for (m in 2:12)
  {
    Lep<-c()
    nam <- paste("Lepage_raw_m", m, sep = "")
    for (n in m:(25-m))
    {
      a<-foreach (i = 1:100000, .combine = 'c') %dopar%
        {
          X<-rnorm(m,0,1)
          Y<-rnorm(n,0,1)
          D <- Lepage_statistic(X,Y)
          D
        }
      Lep <- cbind(Lep,a)
    }
      colnames(Lep)<- paste("n=",c(m:(25-m)))

      assign(nam, Lep)
    
  }
  # Name columns
  ####################
  colnames(Lepage_raw_m10)<- paste("n=",c(10:(25-10)))
  colnames(Lepage_raw_m11)<- paste("n=",c(11:(25-11)))
  colnames(Lepage_raw_m12)<- paste("n=",c(12:(25-12)))
  colnames(Lepage_raw_m2)<- paste("n=",c(2:(25-2)))
  colnames(Lepage_raw_m3)<- paste("n=",c(3:(25-3)))
  colnames(Lepage_raw_m4)<- paste("n=",c(4:(25-4)))
  colnames(Lepage_raw_m5)<- paste("n=",c(5:(25-5)))
  colnames(Lepage_raw_m6)<- paste("n=",c(6:(25-6)))
  colnames(Lepage_raw_m7)<- paste("n=",c(7:(25-7)))
  colnames(Lepage_raw_m8)<- paste("n=",c(8:(25-8)))
  colnames(Lepage_raw_m9)<- paste("n=",c(9:(25-9)))
  ####################
  # Lepage upper tail prob tables
  
  ####################################
  
  final <- matrix(nrow = 2*132, ncol = 4)
  count <-0
  
  for (m in 2:12)
  {
    for(n in m:(25-m))
    {
      count <- count + 1
      final[2*(count)-1, 1] <- paste("m =",m)
      final[2*count, 1] <- paste("m =",m)
      
      final[2*count-1, 2] <- paste("n =",n)
      final[2*count, 2] <- paste("n =",n)
    }
  }
  
  critical_approx <- c(0.05, 0.1)
  
  count <-0
  
  for (index in 2:12)
  {
    l <- dim(get( paste("Lepage_raw_cauchy_m",index,sep = "") )  )[2]
    A <-get( paste("Lepage_raw_m",index,sep = "") )
    tmp <- matrix(nrow = 2*l, ncol = 2 )
    tmp_name <- c()
    for (j in index:(25-index) ) {
      
      
      tmp_name <- c(tmp_name, as.character(rep(j,2)) )
      
      
    }
    
    nam <-paste("Lepage_upper_tail_prob_m",index,sep = "")
   
    # row.names( tmp ) <- paste(" n = ", tmp_name, sep = "  ")
     # colnames(tmp) <-paste("n=",c(3:min(index,25-index)))
     # colnames(final) <- c("m","n","X", "P(D \u2265 X)")
     
    for (c in 1:l)
    {
      vect <- A[,c]
      iter <- 1
      closest_points <- matrix(rep(10,4), byrow = TRUE, nrow = 2, ncol = 2)
      all_points <- matrix( nrow = length(unique(vect)), ncol = 2 )
      for (crit_val in unique(vect)   )
      {
        per <-length( vect[vect>=crit_val])/length(vect)
        # tmp[iter, c] <- per
        all_points[iter,1] <- crit_val
        all_points[iter,2] <- per
        iter<- iter + 1
        
      }
      tmp_points <- all_points[,1]
      tmp_prob <- all_points[,2]
      
      index_5_bigger  <-match(min(tmp_prob[tmp_prob>=0.05]   ),tmp_prob)
      # index_5_smaller <- match(max(tmp_prob[tmp_prob<=0.05]   ),tmp_prob)
      
      index_10_bigger  <-match(min(tmp_prob[tmp_prob>=0.1]   ),tmp_prob)
      # index_10_smaller <- match(max(tmp_prob[tmp_prob<=0.1]   ),tmp_prob)
      
      
      # closest_points[1,]<- c( tmp_points[index_5_smaller], tmp_prob[index_5_smaller] )
      closest_points[1,]<- c( tmp_points[index_5_bigger], tmp_prob[index_5_bigger] )
      # closest_points[3,]<- c( tmp_points[index_10_smaller], tmp_prob[index_10_smaller] )
      closest_points[2,]<- c( tmp_points[index_10_bigger], tmp_prob[index_10_bigger] )
      # tmp[((c-1)*4+1):(c*4)  ,1:2] <-closest_points
      # tmp[c  ,1:2] <-closest_points
      count <- count+1
      final[ count:(count+1),3:4] <- closest_points
      count <- count + 1
      
    }
    # assign(nam,  tmp )
     # final[]
  }
  colnames(final) <- c("m","n","X", "P(D \u2265 X)")
  
  save( file="Lepage_upper_probablity_tables.RData",
        Lepage_upper_tail_prob_m2,
        Lepage_upper_tail_prob_m3,Lepage_upper_tail_prob_m4, Lepage_upper_tail_prob_m5,
        Lepage_upper_tail_prob_m6, Lepage_upper_tail_prob_m7, Lepage_upper_tail_prob_m8,
        Lepage_upper_tail_prob_m9,Lepage_upper_tail_prob_m10,Lepage_upper_tail_prob_m11,
        Lepage_upper_tail_prob_m12)
  
  
  
  save( file="Lepage_original_raw_data.RData",
        Lepage_raw_m2,
        Lepage_raw_m3,Lepage_raw_m4, Lepage_raw_m5,
        Lepage_raw_m6, Lepage_raw_m7, Lepage_raw_m8,
        Lepage_raw_m9,Lepage_raw_m10,Lepage_raw_m11,
        Lepage_raw_m12)
  
  
  
  
  ####################################
  
  
 

  
  
  # Lehman normal data upper prob table
  #####################################
  
  final_Lehman_upper_prob_table <- matrix(nrow = 4*132, ncol = 4)
  count <-0

  
  for (m in 2:12)
  {

    for(n in m:(25-m))
    {
      count <- count + 1
      final_Lehman_upper_prob_table[count, 1] <- paste("m =",m)
      final_Lehman_upper_prob_table[count, 2] <- paste("n =",n)
      
      count <- count + 1
      final_Lehman_upper_prob_table[count, 1] <- paste("m =",m)
      final_Lehman_upper_prob_table[count, 2] <- paste("n =",n)
      
      
      count <- count + 1
      final_Lehman_upper_prob_table[count, 1] <- paste("m =",m)
      final_Lehman_upper_prob_table[count, 2] <- paste("n =",n)
      
      
      count <- count + 1
      final_Lehman_upper_prob_table[count, 1] <- paste("m =",m)
      final_Lehman_upper_prob_table[count, 2] <- paste("n =",n)
      
      
   
    }
  }
  
  critical_approx <- c(0.05, 0.1)
  
  count <-0
  
  for (index in 2:12)
  {
    l <- dim(get( paste("Lehman_raw_cauchy_m",index,sep = "") )  )[2]
    A <-get( paste("Lehman_raw_cauchy_m",index,sep = "") )
    tmp <- matrix(nrow = 4*l, ncol = 2 )
    tmp_name <- c()
    for (j in index:(25-index) ) {
      
      
      tmp_name <- c(tmp_name, as.character(rep(j,2)) )
      
      
    }
    
    nam <-paste("Lehman_upper_tail_prob_m",index,sep = "")
    
    # row.names( tmp ) <- paste(" n = ", tmp_name, sep = "  ")
    # colnames(tmp) <-paste("n=",c(3:min(index,25-index)))
    # colnames(final) <- c("m","n","X", "P(D \u2265 X)")
    
    for (c in 1:l)
    {
      vect <- A[,c]
      iter <- 1
      closest_points <- matrix(rep(10,8), byrow = TRUE, nrow = 4, ncol = 2)
      all_points <- matrix( nrow = length(unique(vect)), ncol = 2 )
      for (crit_val in unique(vect)   )
      {
        per <-length( vect[vect>=crit_val])/length(vect)
        # tmp[iter, c] <- per
        all_points[iter,1] <- crit_val
        all_points[iter,2] <- per
        iter<- iter + 1
        
      }
      tmp_points <- all_points[,1]
      tmp_prob <- all_points[,2]
      
      index_5_bigger  <-match(min(tmp_prob[tmp_prob>=0.025]   ),tmp_prob)
      index_5_smaller <- match(max(tmp_prob[tmp_prob<=0.025]   ),tmp_prob)
      
      index_10_bigger  <-match(min(tmp_prob[tmp_prob>=0.975]   ),tmp_prob)
      index_10_smaller <- match( max(tmp_prob[tmp_prob<=0.975]   ),tmp_prob)
      
      
      closest_points[1,]<- c( tmp_points[index_5_smaller], tmp_prob[index_5_smaller] )
      closest_points[2,]<- c( tmp_points[index_5_bigger], tmp_prob[index_5_bigger] )
      closest_points[3,]<- c( tmp_points[index_10_smaller], tmp_prob[index_10_smaller] )
      closest_points[4,]<- c( tmp_points[index_10_bigger], tmp_prob[index_10_bigger] )
      # tmp[((c-1)*4+1):(c*4)  ,1:2] <-closest_points
      # tmp[c  ,1:2] <-closest_points
      count <- count+1
      final_Lehman_upper_prob_table[ count:(count+3),3:4] <- closest_points
      count <- count +3
      
    }
    # assign(nam,  tmp )
    # final[]
  }
  colnames(final_Lehman_upper_prob_table) <- c("m","n","x", "P(V \u2265 x)")
  final_Lehman_upper_prob_table <- na.omit(final_Lehman_upper_prob_table)
  write.csv2(final_Lehman_upper_prob_table, file =  "final_Lehman_upper_prob_table_cauchy.csv")
  
  #####################################
  