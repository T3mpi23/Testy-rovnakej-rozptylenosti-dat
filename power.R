
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






      replic <-10000
      
          ay_Fisher <- c()
          ay_Mood <-c()
          ay_Ansari<-c()
          ay_Lehman <-c()
          ay_Miller <-c()
          ay_Lehman_pivot <-c()
          ay_Lehman_perm <-c()
          ay_Lepage <- c()
          
          m = 30#size of X data
          n = 30# size of Y data   Ansari n <= m   (m + n)<= 20
         
       
           s <- m + n
         # kritiké hodnoty Ansariho-Bradleyho test. štatistiky
           c_AB_upper <- 68  
          c_AB_lower <- 42 - 1
          
           
        
           d_lepage_critical <- 5.718615



           student_criticals <- qt(0.975,df=m+n-2)
        
           critical <- qnorm(0.975,mean=0,sd=1)
        chi_critical <- qchisq(0.95, 2)
        
          mx <-0
        my <-0
        
        sdx_max <-2
        sdy<-1
        tick<-0.1
        osX<-seq(0.40, sdx_max,by=tick)
        # poèet kôl bootstrapu (resp. permutovania)
          boot_rep <-200
         
  
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
                
                # Normal
 

                X <- rnorm(m,mean = mx, sd = sqrt(sdx))
                Y <- rnorm(n, mean = my, sd = sqrt(1))
               # cauchy
                # X <- rcauchy(m, location = mx,scale = sdx )
                # Y <- rcauchy(n, location = my,scale = 1)
              # student
                  # X <- rt(m, df= 5)*(sqrt(sdx))
                # Y <- rt(n, df = 5)
               
                #FISHER TEST
                #############################
                if (var.test(X,Y, conf.level = 0.95)$p.value <= 0.05)
                {
                count_Fisher <- count_Fisher + 1
                }
                # ##############################
              
                #Moood
                #############################

                W <- Mood_statistic(X,Y)
                
                E_W <-n*(s+1)*(s-1)/12
                D_W <- m*n*(s+1)*(s+2)*(s-2)/180
                if ( ( abs(W-E_W)/sqrt(D_W) >= critical ) )
                  {
                   count_Mood <- count_Mood + 1
                  
                }
                        
        
                ##############################
                
                # ANSARI-Bradley blok
                #########################################
                C <- Ansari_statistic(X,Y)
               
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

                  if ( abs( ( C-E_C )/sqrt(D_C) ) >= critical )
                  {
                   
                   count_Ansari <- count_Ansari + 1
                  }
                ##########################################
                # LEHMAN
                ######################################################
             
                L<-Lehman_statistic(X,Y)
                E_L <- m*(m-1)*n*(n-1)/8
                D_L_perm <- var (Permutaion_test(X,Y,boot_rep))

                if(  (  abs ((L-E_L)/sqrt(D_L_perm))  >= critical )  )
                {
                 count_Lehman_perm <- count_Lehman_perm+ 1 
                }
  
               
                D_L_boot <- var(bootstrap(X ,Y ,boot_rep))
                if(  (  abs ((L-E_L)/sqrt(D_L_boot))  >= critical )  )
                {
                  count_Lehman_asym = count_Lehman_asym + 1
                }

                # Pivotal Lehman
                L_boot <- bootstrap(X,Y,boot_rep)
              
                if( L<= quantile(L_boot, probs = 0.025) | L>= quantile(L_boot, probs = 0.975))
                {
                count_Lehman_pivot <- count_Lehman_pivot + 1
                }

               
                ######################################################
                
                #Miller test
                ###################
                
                Q <- Miller_statistic(X,Y)
                if ( abs(Q)>=critical )
                  {
                 
                   count_Miller <- count_Miller + 1
                  
                }

                ###################
                  #  LEPAGE
                # ###################
                D <- Lepage_statistic(Y,X)
                
                  if (D >= chi_critical){
                    count_Lepage <- count_Lepage + 1

                
               } #end of for
               
             
          
             ay_Fisher[length(ay_Fisher)] <- count_Fisher/replic
             ay_Mood[length(ay_Mood)] <- count_Mood/replic
             ay_Ansari[length(ay_Ansari)] <- count_Ansari/replic
            ay_Lehman[length(ay_Lehman)] <- count_Lehman_asym/replic
            ay_Lehman_perm[length(ay_Lehman_perm)] <- count_Lehman_perm / replic
            ay_Lehman_pivot[length(ay_Lehman_pivot)] <- count_Lehman_pivot / replic
            ay_Miller[length(ay_Miller)] <- count_Miller/replic
            ay_Lepage[length(ay_Lepage)] <- count_Lepage/replic
              
              }
              
                }
            d2 <- Sys.time()
            d2 - d1
        
            
             
        

   # Lepage tables
  #########################
  
   Lep<-c()
   for (m in 2:12)
   {
     Lep<-c()
     nam <- paste("Lep", m, sep = "")
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
      
      
      closest_points[1,]<- c( tmp_points[index_5_bigger], tmp_prob[index_5_bigger] )
      closest_points[2,]<- c( tmp_points[index_10_bigger], tmp_prob[index_10_bigger] )
      count <- count+1
      final[ count:(count+1),3:4] <- closest_points
      count <- count + 1
      
    }
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
   
      count <- count+1
      final_Lehman_upper_prob_table[ count:(count+3),3:4] <- closest_points
      count <- count +3
      
    }
    
  }
  colnames(final_Lehman_upper_prob_table) <- c("m","n","x", "P(V \u2265 x)")
  final_Lehman_upper_prob_table <- na.omit(final_Lehman_upper_prob_table)
  write.csv2(final_Lehman_upper_prob_table, file =  "final_Lehman_upper_prob_table_cauchy.csv")
  
  #####################################
  