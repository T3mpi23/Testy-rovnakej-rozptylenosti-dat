
# Moodova testova štatistika
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

# Ansariho-Bradleyho testova štatistika
Ansari_statistic <- function(X,Y)
{
  rank_vector <-c()
  even <- 0
  m <- length(X)
  n <- length(Y)
  s <- m+n
  C <-0
   rank_vector <- rank(c(X,Y))  #spojena sad dat
  
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

#Lepageova testova štatistika
Lepage_statistic <- function(X,Y )
{
  m <- length(X)
  n <- length(Y)
  
  #Wilcoxonova testova štatistika
  W = sum(rank(c(Y,X))[1:n] ) 
  E_W <-n*(m+n+1)/2  #stedna hodnota E(W)
  D_W <- m*n*(m+n+1)/12 #rozptyl D(W)

  s <- m+n
   
  #Stredna hodnota a rozptyl Ansariho-Bradleyho test. štat.
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

#Millerova testova štatistika
Miller_statistic <- function(X,Y )
{
  m <- length(X)
  n <- length(Y)
  
  #pomocne matice, kde kazdy riadok je sada dat X (resp. Y)
  #Odstranenim diagonaly ziskame sady dat v ktorych postupne 
  #vynechavame prvky (Jacknife)
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
#Lehmanova testova štatistika 
Lehman_statistic <-function(X,Y){
  
   m<- length(X)
   n<- length(Y)
   
   # Kroneckerov produkt pre urychlenie, urobi všetky rozdiely
    X_diffs <-  matrix ( abs(kronecker(X,X,FUN ="-")), 
                         byrow = TRUE, nrow = m)
    X_diffs<- X_diffs[ upper.tri(X_diffs)]
    
    Y_diffs  <-  matrix ( (abs(kronecker(Y,Y,FUN ="-")  )), 
                          byrow = TRUE, nrow = n)
    Y_diffs <- Y_diffs[upper.tri(Y_diffs)]
    # rozdiely rozdielov, berieme len pozitivne
   M <- kronecker(X_diffs,Y_diffs, FUN = "-" )[kronecker(X_diffs,
                                                         Y_diffs, FUN = "-")>0] 
   #Lehmannova test. štat.
   L <- length(   M   )

   return(L)
}
# Permutacny pristup odhadu disperzie Lehmannovej test. štat.
# B je pocet kôl
# vystup: vektor Lehmanovych test. štat.
# odhad disperzie je potom rozptyl vystupu
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
#B pocet kôl
# vystup: vektor Lehmanovych test. štat.
bootstrap  <- function(X,Y,B) 
  {
  L_boot <-c()
  X_rand <-c()
  Y_rand <-c()
  m <- length(X)
  n <- length(Y)
  for (b in  1:B)
  {
    # tvorba novych sad dat
    rand_index <- sample(1:( m+n), length(X)+length(Y), replace =  TRUE)
      X_rand <- c(X,Y)[rand_index[1:m]]
      Y_rand <- c(X,Y)[rand_index[(m+1):(n+m)]]
    L_boot[length(L_boot)+1] <- Lehman_statistic(X_rand, Y_rand)
  }
  return( (L_boot) )
}
# modifikovany bootstrap kde medzi sebou nemiešame prvky z rôznych sad
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
# m,n velkosti sád dát (m <= n)
# replic - pocet simulacii pre jednotlive pomery rozptylov
# mx, my - stredne hodnoty rozdeleni z ktorych su generovane data
# B pocet kol bootstrapu, resp. permutovania
# dist - premenna ktora urcuje z akeho rozdelenia generujeme data:
        # dist = -1 : normalne rozdelenie
        # dist = 0  : Cauchyho rozdelenie
        # dist > 0  : Studentovo t rozdelenie s dist stupnami volnosti
Power_function_estimate <- function(m,n, replic, mx,my, B, dist)
{
          ay_Fisher <- c()
          ay_Mood <-c()
          ay_Ansari<-c()
          ay_Lehman <-c()
          ay_Miller <-c()
          ay_Lehman_pivot <-c()
          ay_Lehman_perm <-c()
          ay_Lepage <- c()
          s <- m + n
          critical <- qnorm(0.975,mean=0,sd=1)
          chi_critical <- qchisq(0.95, 2)
         
          osX<-seq(0.40, 3,by=0,1)
        # poèet kôl bootstrapu (resp. permutovania)
          boot_rep <-B
             for (sdx in osX)
              {
              count_Fisher <- 0
              count_Mood <- 0
              count_Ansari <- 0
              count_Lehman_asym <- 0
              count_Lehman_perm <- 0
              count_Lehman_pivot <- 0 
              count_Miller <- 0
              count_Lepage <- 0
              for( i in 1:replic)
                  {
                # Normal data generating
                if (dist == 0)
                {
                   X <- rcauchy(m, location = mx,scale = sdx )
                   Y <- rcauchy(n, location = my,scale = 1)
                }
                if (dist == -1)
                {
                   X <- rnorm(m,mean = mx, sd = sqrt(sdx))
                   Y <- rnorm(n, mean = my, sd = sqrt(1))
                }
                
                if (dist > 0)
                {
                  X <- rt(m, df= dist)*(sqrt(sdx)) + mx
                  Y <- rt(n, df = dist) + my
                }
             
                #FISHER TEST
                if (var.test(X,Y, conf.level = 0.95)$p.value <= 0.05)  {
                  count_Fisher <- count_Fisher + 1   }
                
            
                #Moood
                W <- Mood_statistic(X,Y)
                E_W <-n*(s+1)*(s-1)/12
                D_W <- m*n*(s+1)*(s+2)*(s-2)/180
                if ( ( abs(W-E_W)/sqrt(D_W) >= critical ) ){
                   count_Mood <- count_Mood + 1   }
              
                # ANSARI-Bradley blok
                C <- Ansari_statistic(X,Y)
               
                if (s %%2==0) {
                  E_C <- n*(s+2)/4
                  D_C <- m*n*(s+2)*(s-2)/(48*(s-1)) }
                if (s %%2==1){
                  E_C <- n*(s+1)^2 / (4*s)
                  D_C <- m*n*(s+1)*(s^2 +3)/(48*(s^2))   }

                  if ( abs( ( C-E_C )/sqrt(D_C) ) >= critical ) {
                   count_Ansari <- count_Ansari + 1    }
               
                 # LEHMAN
                L<-Lehman_statistic(X,Y)
                E_L <- m*(m-1)*n*(n-1)/8
                D_L_perm <- var (Permutaion_test(X,Y,boot_rep))

                if(  (  abs ((L-E_L)/sqrt(D_L_perm))  >= critical )  ) {
                 count_Lehman_perm <- count_Lehman_perm+ 1  }
  
                D_L_boot <- var(bootstrap(X ,Y ,boot_rep))
                if(  (  abs ((L-E_L)/sqrt(D_L_boot))  >= critical )  ) {
                  count_Lehman_asym = count_Lehman_asym + 1 }

                # Pivotal Lehman
                L_boot <- bootstrap(X,Y,boot_rep)
              
                if( L<= quantile(L_boot, probs = 0.025) | L>= quantile(L_boot, 
                                                              probs = 0.975)){
                count_Lehman_pivot <- count_Lehman_pivot + 1  }

                #Miller test
                Q <- Miller_statistic(X,Y)
                if ( abs(Q)>=critical ) {count_Miller <- count_Miller + 1 }

                  #  LEPAGE
                D <- Lepage_statistic(Y,X)
                
                  if (D >= chi_critical){count_Lepage <- count_Lepage + 1 } 
                #end of for
              }
              
             ay_Fisher[length(ay_Fisher)] <- count_Fisher/replic
             ay_Mood[length(ay_Mood)] <- count_Mood/replic
             ay_Ansari[length(ay_Ansari)] <- count_Ansari/replic
            ay_Lehman[length(ay_Lehman)] <- count_Lehman_asym/replic
            ay_Lehman_perm[length(ay_Lehman_perm)] <- count_Lehman_perm / replic
            ay_Lehman_pivot[length(ay_Lehman_pivot)] <- count_Lehman_pivot / replic
            ay_Miller[length(ay_Miller)] <- count_Miller/replic
            ay_Lepage[length(ay_Lepage)] <- count_Lepage/replic
              }
              
              output <- matrix(ncol = length(osX), nrow = 9, byrow = TRUE,osX,
                               ay_Fisher, ay_Mood,ay_Ansari,ay_Lehman,
                               ay_Lehman_pivot,ay_Lehman_perm,ay_Miller,
                               ay_Lepage  )
        rownames(output)  <- c( "pomer rozptylov" , 
                               "Fisher", "Mood","Ansari-Bradley", 
                               "Lehman- bootstrap", "Lehman-pivot ", 
                               "Lehman-permutation", "Miller", "Lepage")
        return(output)
             }

     
            # vytvorenie clusterov pre paralelny program,
            # znacne urychli vypocet
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
            n.cores <- parallel::detectCores() - 5
            
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
            

  
Lepage_tables <- function()
{  
  ####### Lepage raw data generating
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
  ###################
  
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
  
  
  
  count <-0
  
  for (index in 2:12)
  {
    l <- dim(get( paste("Lepage_raw_m",index,sep = "") )  )[2]
    A <-get( paste("Lepage_raw_m",index,sep = "") )
    tmp <- matrix(nrow = 2*l, ncol = 2 )
    tmp_name <- c()
    for (j in index:(25-index) ) {
      
      
      tmp_name <- c(tmp_name, as.character(rep(j,2)) )
      
      
    }
    
    nam <-paste("Lepage_upper_tail_prob_m",index,sep = "")
   
    
     
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

      index_10_bigger  <-match(min(tmp_prob[tmp_prob>=0.1]   ),tmp_prob)

      
      closest_points[1,]<- c( tmp_points[index_5_bigger], tmp_prob[index_5_bigger] )
      closest_points[2,]<- c( tmp_points[index_10_bigger], tmp_prob[index_10_bigger] )
      count <- count+1
      final[ count:(count+1),3:4] <- closest_points
      count <- count + 1
      
    }
  }
  colnames(final) <- c("m","n","X", "P(D \u2265 X)")

  ####################################
  return(final)
}
 

Lehman_cauchy_tables <- function()
{
  

# Lehman raw data cauchy
  set.seed(1)
  for (m in 2:12)
  {
    Lep<-c()
    nam <- paste("Lehman_raw_cauchy_m", m, sep = "")
    for (n in m:(25-m))
    {
      a<-foreach (i = 1:100000, .combine = 'c') %dopar%
        {
          X<-rcauchy(m,0,1)
          Y<-rcauchy(n,0,1)
          V <- Lehman_statistic(X,Y)
          V
        }
      Lep <- cbind(Lep,a)
    }
    colnames(Lep)<- paste("n=",c(m:(25-m)))
    
    assign(nam, Lep)
    
  }
  
  
  # Lehman cauchy data upper prob table
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

    for (c in 1:l)
    {
      vect <- A[,c]
      iter <- 1
      closest_points <- matrix(rep(10,8), byrow = TRUE, nrow = 4, ncol = 2)
      all_points <- matrix( nrow = length(unique(vect)), ncol = 2 )
      for (crit_val in unique(vect)   )
      {
        per <-length( vect[vect>=crit_val])/length(vect)
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
      
      count <- count+1
      final_Lehman_upper_prob_table[ count:(count+3),3:4] <- closest_points
      count <- count +3
      
    }
    
  }
  
  colnames(final_Lehman_upper_prob_table) <- c("m","n","x", "P(V \u2265 x)")
  final_Lehman_upper_prob_table_cauchy <- na.omit(final_Lehman_upper_prob_table)
  return(final_Lehman_upper_prob_table_cauchy)
}

Lehman_normal_tables <- function()
{
  

  # Lehman raw data normal
  set.seed(1)
  for (m in 2:12)
  {
    Lep<-c()
    nam <- paste("Lehman_raw_m", m, sep = "")
    for (n in m:(25-m))
    {
      a<-foreach (i = 1:100000, .combine = 'c') %dopar%
        {
          X<-rnorm(m,0,1)
          Y<-rnorm(n,0,1)
          V <- Lehman_statistic(X,Y)
          V
        }
      Lep <- cbind(Lep,a)
    }
    colnames(Lep)<- paste("n=",c(m:(25-m)))
    
    assign(nam, Lep)
    
  }
  
  
  
  
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
  count <-0
  
  for (index in 2:12)
  {
    l <- dim(get( paste("Lehman_raw_m",index,sep = "") )  )[2]
    A <-get( paste("Lehman_raw_m",index,sep = "") )
    tmp <- matrix(nrow = 4*l, ncol = 2 )
    tmp_name <- c()
    for (j in index:(25-index) ) {
      tmp_name <- c(tmp_name, as.character(rep(j,2)) )
    }
    
    nam <-paste("Lehman_upper_tail_prob_m",index,sep = "")
    
    for (c in 1:l)
    {
      vect <- A[,c]
      iter <- 1
      closest_points <- matrix(rep(10,8), byrow = TRUE, nrow = 4, ncol = 2)
      all_points <- matrix( nrow = length(unique(vect)), ncol = 2 )
      for (crit_val in unique(vect)   )
      {
        per <-length( vect[vect>=crit_val])/length(vect)
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
      
      count <- count+1
      final_Lehman_upper_prob_table[ count:(count+3),3:4] <- closest_points
      count <- count +3
      
    }
    
  }
  
  colnames(final_Lehman_upper_prob_table) <- c("m","n","x", "P(V \u2265 x)")
  final_Lehman_upper_prob_table <- na.omit(final_Lehman_upper_prob_table)
return(final_Lehman_upper_prob_table)
}
  #####################################

