args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Provide a b", call.=FALSE)
}

#Functions

hiveq <- function(d, BirthRate, b, Sf) {
  write.csv(d,
            file = paste0(deparse(BirthRate),".",deparse(b),".",deparse(Sf),".","Figure6VeryLongBenefit.csv"),
            row.names = FALSE)
}

fib <- function(n){
  fib <- vector()
  fib[1] <- 1
  fib[2] <- 1
  for(i in 3:n){
    fib[i] <- fib[i-1] + fib[i-2]
  }
  return(fib[n])
}

printsummary <- function(Gen, means, sds) {
  cat("At Gen", Gen," Mean starting male prop is:",means[1]," Mean finishing male prop is:",means[2]," Mean helping tendency is:",means[3]," Standard Dev starting male prop is:",sds[1]," Standard Dev finishing male prop is:",sds[2]," Standard Dev Helping Tendency is:",sds[3])
} #This Function Prints the Mean Spring Female Ratio, Summer Female Ratio, and Helper Ratio at a Given Generation

mutate <- function(m){
  if (Gen > Lag | m!=0){
    if (runif(1) < mu){ #If below the mutation rate per allele
      m = m + rnorm(1, mean=0, sd=sdmu) #Change the value by that much
    }
    if(m<0){ #Can't have values less than 0 or greater than 1
      m = 0
    }
    if(m>1){
      m = 1
    }
  }
  return(m)
} #This Function Mutates a value with prob. mu by a stdev of sdmu

mutate2 <- function(m){
  if (Gen > Lag | m!=0){
    if (runif(1) < mu){ #If below the mutation rate per allele
      m = m + rnorm(1, mean=0, sd=sdmu2) #Change the value by that much
    }
    if(m< -0.01){ #Can't have values less than -0.01 or greater than 0.01
      m = -0.01
    }
    if(m> 0.01){
      m = 0.01
    }
  }
  return(m)
} #This Function Mutates a value with prob. mu by a stdev of sdmu (I put this twice so I can run the code twice)

Stats <- function(F1){
  N <- as.numeric(length(F1))
  Totals = c(0, 0, 0, 0)
  ss = c(0, 0, 0, 0)
  means = c(0, 0, 0, 0)
  sd = c(0, 0, 0, 0)
  for(i in 1:N){ 
    z1 <- 0.5 * (F1[[i]][["Female"]][["Starting Male"]][1] + F1[[i]][["Female"]][["Starting Male"]][2]) #Average of male ratio (chrom 1) with male ratio (chrom 2) at start
    Totals[1] <- Totals[1] + z1
    ss[1] <- ss[1] + (z1*z1)
    z2 <- 0.5 * (F1[[i]][["Female"]][["Finishing Male"]][1] + F1[[i]][["Female"]][["Finishing Male"]][2]) #Average of male ratio (chrom 1) with malle ratio (chrom 2) at end
    Totals[2] <- Totals[2] + z2
    ss[2] <- ss[2] + (z2*z2)
    h1 <-  0.5 * (F1[[i]][["Female"]][["Starting Hel"]][1] + F1[[i]][["Female"]][["Starting Hel"]][2]) #Sum of helpers (chrom) with helpers (chrom 2)
    Totals[3] <- Totals[3] + h1
    ss[3] <- ss[3] + (h1*h1)
    h2 <-  0.5 * (F1[[i]][["Female"]][["Finishing Hel"]][1] + F1[[i]][["Female"]][["Finishing Hel"]][2]) #Sum of helper finish (chrom 1) with helper finish (chrom 2)
    Totals[4] <- Totals[4] + h2
    ss[4] <- ss[4] + (h2*h2)
  }
  means[1] <- Totals[1]/N
  sd[1] <- sqrt(ss[1]/(N*(means[1]*means[1])))
  means[2] <- Totals[2]/N
  sd[2] <- sqrt(ss[2]/(N*(means[2]*means[2])))
  means[3] <- Totals[3]/N
  sd[3] <- sqrt(ss[3]/(N*(means[3]*means[3])))
  means[4] <- Totals[4]/N
  sd[4] <- sqrt(ss[4]/(N*(means[4]*means[4])))
  Data <- c(means, sd)
  return(Data)
} #This Function calculates the Stats for a set of nests

NextGeneration <- function(F1){
  Son <- Male1
  Daughter <- Female1
  NF2 <- 1
  NF3 <- as.numeric(length(F1)) + 1
  NM1 <- 1
  NM2 <- 1
  F3[1:(NF3-1)] <- F1 
  for(j in 1:DIY){
    if(NF2 > 1){
      for(i in 1:(NF2-1)){
        if(F2[[i]][["Age"]] == RepAge + TimeToNest){
          if(NM2 > 5){
            F3[[NF3]][["Female"]] <- F2[[i]][["Female"]]
            F3[[NF3]][["Male"]] <- M2[[sample(1:(NM2-2), 1)]]
            NF3 = NF3 + 1
            F2[[i]] = F2[[NF2-1]]
            NF2 = NF2 - 1
          }
          else{
            F2[[i]][["Age"]] = F2[[i]][["Age"]] - 2 #If no males are mature, hold off
          }
        }   #Fertilize Mature Reproductive Females
        else{
          if(runif(1) < Sfe){
            F2[[i]] = F2[[NF2-1]]
            NF2 = NF2 - 1
            F2[[i]][["Age"]] = F2[[i]][["Age"]] + 1
          }
          else{
            F2[[i]][["Age"]] = F2[[i]][["Age"]] + 1
          }
        }
      }
    }
    for(i in 1:(NF3-1)){
      Kids = rpois(1, F3[[i]][["SBR"]])
      if(Kids > 0){
        for(p in 1:Kids){
          z1 <- 0.5 * (F3[[i]][["Female"]][["Starting Male"]][1] + F3[[i]][["Female"]][["Starting Male"]][2]) #Average of male (chrom 1) with male (chrom 2) starting ratios
          z2 <- 0.5 * (F3[[i]][["Female"]][["Finishing Male"]][1] + F3[[i]][["Female"]][["Finishing Male"]][2]) #Average of male (chrom 1) with male (chrom 2) finishing ratios
          z <- (as.numeric(j))*(z2) + z1
          if(z<0){ #Can't have values less than 0 or greater than 1
            z <- 0
          }
          if(z>1){
            z <- 1
          }
          if(runif(1) < z){ #If Male (remember z is male rate)
            if(LD == 0){
              SonMSta <- F3[[i]][["Female"]][["Starting Male"]][sample(1:2, 1)] #Pick a Random Chromosome from Female (Haplodiploid!)
              SonMFin <- F3[[i]][["Female"]][["Finishing Male"]][sample(1:2, 1)] #Same with the Finishing Prop Rate
              SonHSta <- F3[[i]][["Female"]][["Starting Hel"]][sample(1:2, 1)] #Same with the Helping Rate
              SonHFin <- F3[[i]][["Female"]][["Finishing Hel"]][sample(1:2, 1)] #Same with the Helping Rate
            }
            else{
              Chr <- sample(1:2, 1)
              SonMSta <- F3[[i]][["Female"]][["Starting Male"]][Chr] #Pick a Random Chromosome from Female (Haplodiploid!)
              SonMFin <- F3[[i]][["Female"]][["Finishing Male"]][Chr] #Same with the Finishing Prop Rate
              SonHSta <- F3[[i]][["Female"]][["Starting Hel"]][Chr] #Same with the Helping Rate
              SonHFin <- F3[[i]][["Female"]][["Finishing Hel"]][Chr] #Same with the Helping Rate
            }
            Son[["Starting Male"]] <- mutate(SonMSta) #MutateTheValuesAndCreateASonWithThem
            Son[["Finishing Male"]] <- mutate2(SonMFin)
            Son[["Starting Hel"]] <- mutate(SonHSta)
            Son[["Finishing Hel"]] <- mutate2(SonHFin)
            Son[["Age"]] <- 0
            M1[[NM1]] <- Son
            NM1 <- NM1 + 1
          }
          else{
            if(LD == 0){
              DaughterMSta <- F3[[i]][["Female"]][["Starting Male"]][sample(1:2, 1)] #Pick a Random Chromosome from Female (Haplodiploid!)
              DaughterMFin <- F3[[i]][["Female"]][["Finishing Male"]][sample(1:2, 1)] #Same with the Finishing Prop Rate
              DaughterHSta <- F3[[i]][["Female"]][["Starting Hel"]][sample(1:2, 1)] #Same with the Helping Rate
              DaughterHFin <- F3[[i]][["Female"]][["Finishing Hel"]][sample(1:2, 1)] #Same with the Helping Rate
            }
            else{
              Chr <- sample(1:2, 1)
              DaughterMSta <- F3[[i]][["Female"]][["Starting Male"]][Chr] #Pick a Random Chromosome from Female (Haplodiploid!)
              DaughterMFin <- F3[[i]][["Female"]][["Finishing Male"]][Chr] #Same with the Finishing Prop Rate
              DaughterHSta <- F3[[i]][["Female"]][["Starting Hel"]][Chr] #Same with the Helping Rate
              DaughterHFin <- F3[[i]][["Female"]][["Finishing Hel"]][Chr] #Same with the Helping Rate
            }
            Daughter[["Starting Male"]][1] <- mutate(DaughterMSta) #MutateTheValuesAndCreateASonWithThem
            Daughter[["Finishing Male"]][1] <- mutate2(DaughterMFin)
            Daughter[["Starting Hel"]][1] <- mutate(DaughterHSta)
            Daughter[["Finishing Hel"]][1] <- mutate2(DaughterHFin)
            DaughterMSta <- F3[[i]][["Male"]][["Starting Male"]]
            DaughterMFin <- F3[[i]][["Male"]][["Finishing Male"]]
            DaughterHSta <- F3[[i]][["Male"]][["Starting Hel"]]
            DaughterHFin <- F3[[i]][["Male"]][["Finishing Hel"]]
            Daughter[["Starting Male"]][2] <- mutate(DaughterMSta) #MutateTheValuesAndCreateASonWithThem
            Daughter[["Finishing Male"]][2] <- mutate2(DaughterMFin)
            Daughter[["Starting Hel"]][2] <- mutate(DaughterHSta)
            Daughter[["Finishing Hel"]][2] <- mutate2(DaughterHFin)
            h1=0.5*(Daughter[["Starting Hel"]][1]+Daughter[["Starting Hel"]][2])
            h2=0.5*(Daughter[["Finishing Hel"]][1]+Daughter[["Finishing Hel"]][2])
            h <- (as.numeric(j))*(h2) + h1
            if(h<0){ #Can't have values less than 0 or greater than 1
              h = 0
            }
            if(h>1){
              h = 1
            }
            if(runif(1) < h){
              for (k in 1:RepAge){
                if (is.na(F3[[i]][["HelAge"]][k])){
                  F3[[i]][["HelAge"]][k] = 0
                  break
                }
              }
            }
            else{
              F2[[NF2]][["Female"]] <- Daughter
              F2[[NF2]][["SBR"]] <- BirthRate
              F2[[NF2]][["HelAge"]] <- rep(NA, RepAge)
              F2[[NF2]][["Age"]] <- 0
              NF2 = NF2 + 1
            }
          }
        }
      }
      if(runif(1) < Sfe & NF3 > 50){
        F3[[i]] = F3[[NF3-1]]
        NF3 = NF3 - 1
        F3[[i]][["Age"]] <- F3[[i]][["Age"]] + 1
        F3[[i]][["HelAge"]] = F3[[i]][["HelAge"]] + 1/RepAge
        X = which(F3[[i]][["HelAge"]] >= 1)
        if(length(X) != 0){
          F3[[i]][["HelAge"]][X] = NA
          F3[[i]][["SBR"]] <- F3[[i]][["SBR"]] + length(X)*b
        }
      }
      else{
        F3[[i]][["Age"]] <- F3[[i]][["Age"]] + 1
        F3[[i]][["HelAge"]] = F3[[i]][["HelAge"]] + 1/RepAge
        X = which(F3[[i]][["HelAge"]] >= 1)
        if(length(X) != 0){
          F3[[i]][["HelAge"]][X] = NA
          F3[[i]][["SBR"]] <- F3[[i]][["SBR"]] + length(X)*b
        }
      }
    }
    if(NM1 > 1){
      for(i in 1:(NM1-1)){
        if (M1[[i]][["Age"]] == RepAge){
          M2[[NM2]] = M1[[i]] 
          NM2 = NM2 + 1
          M1[[i]] = M1[[NM1-1]]
          NM1 = NM1 - 1
        }
        else{
          if(runif(1) < Sme){
            M1[[i]] = M1[[NM1-1]]
            NM1 = NM1 - 1
            M1[[i]][["Age"]] = M1[[i]][["Age"]] + 1
          }
          else{
            M1[[i]][["Age"]] = M1[[i]][["Age"]] + 1
          }
        }
      }
    }
    if(NM2 > 1){
      for(i in 1:(NM2-1)){      
        if(runif(1) < Sme){
          M2[[i]] = M2[[NM2-1]]
          NM2 = NM2 - 1
        }
      }
    }
  }
  for(i in 1:N){
    RF =  as.numeric(sample(1:(NF3-1), 1))
    F1[[i]] = F3[[RF]]
    F1[[i]][["SBR"]] <- BirthRate
    F1[[i]][["HelAge"]] <- rep(NA, RepAge)
    NF3 = NF3-1
    F3[[RF]] <- F3[[NF3]]
  }
  return(F1)
}

Perturb <- function(F1, SM, FM, SH, FH, Percentage){
  Intruders <- floor(length(F1)*Percentage)
  for(i in 1:Intruders){
    F1[[i]][["Female"]][["Starting Male"]][1] <- SM
    F1[[i]][["Female"]][["Starting Male"]][2] <- SM
    F1[[i]][["Male"]][["Starting Male"]] <- SM
    F1[[i]][["Female"]][["Finishing Male"]][1] <- FM
    F1[[i]][["Female"]][["Finishing Male"]][2] <- FM
    F1[[i]][["Male"]][["Finishing Male"]] <- FM
    F1[[i]][["Female"]][["Starting Hel"]][1] <- SH
    F1[[i]][["Female"]][["Starting Hel"]][2] <- SH
    F1[[i]][["Male"]][["Starting Hel"]] <- SH
    F1[[i]][["Female"]][["Finishing Hel"]][1] <- FH
    F1[[i]][["Female"]][["Finishing Hel"]][2] <- FH
    F1[[i]][["Male"]][["Finishing Hel"]] <- FH
  }
  return(F1)
}

#Set Variables
Reps = as.numeric(args[1])
set.seed(Reps)
b = 0.1 #Increase In Summer Growth Rate from Helpers
BirthRate = 0.05
Sf = 0.1 #Proportion of Foundresses that Die before end of summer
Sm = 0.4 #Mort Rate of Males
DIY <- as.numeric(args[2]) #Length from diapause > diapause
Sfe = Sf / 180 #Per Day Mort Rate
Sme = Sm / 180 #Per Day Mort Rate
mu = 0.05 #Mut Rate Per Allele Per Generation
sdmu = 0.05 #St. Dev. of Mutation Size
sdmu2 = 0.0005
Max = 1000 #MaximumNumberInABrood
N = 50; #Nests
RepAge <- 27
TimeToNest <- 10
NumGen = 2000 #Number Gens
Lag = 0 #Time before helper starts evolving
Skip = 1 #Output Every 10th
MalS01 = 0.5
MalE01 = 0
HelS01 = 0.5
HelE01 = 0
MalS02 = 0.5
MalE02 = 0
HelS02 = 0.5
HelE02 = 0
zFMS1 = c(MalS01, MalS01)
zFMF1 = c(MalE01, MalE01)
zMMS1 = MalS01
zMMF1 = MalE01
hFS1 = c(HelS01, HelS01)
hFF1 = c(HelE01, HelE01)
hMS1 = HelS01
hMF1 = HelE01
zFMS2 = c(MalS02, MalS02)
zFMF2 = c(MalE02, MalE02)
zMMS2 = MalS02
zMMF2 = MalE02
hFS2 = c(HelS02, HelS02)
hFF2 = c(HelE02, HelE02)
hMS2 = HelS02
hMF2 = HelE02
SolSocRatio = 1
SBR = BirthRate
PertNum = 100000000
LD = 0
ChangeInClimate = 12077021
DaysInYearSd = 20222202

#Create Data Frame
Female1 = list()
Female1[[1]] <- zFMS1
Female1[[2]] <- zFMF1
Female1[[3]] <- hFS1
Female1[[4]] <- hFF1
names(Female1) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel")
Male1 = list()
Male1[[1]] <- zMMS1
Male1[[2]] <- zMMF1
Male1[[3]] <- hMS1
Male1[[4]] <- hMF1
Male1[[5]] <- 0
Nest1 = list()
names(Male1) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel", "Age")
Nest1[[1]] = Female1
Nest1[[2]] = Male1
Nest1[[3]] = SBR
Nest1[[4]] = 0
Nest1[[5]] =  rep(NA, RepAge)
names(Nest1) <- c("Female", "Male", "SBR", "Age", "HelAge")
Female2 = list()
Female2[[1]] <- zFMS2
Female2[[2]] <- zFMF2
Female2[[3]] <- hFS2
Female2[[4]] <- hFF2
names(Female2) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel")
Male2 = list()
Male2[[1]] <- zMMS2
Male2[[2]] <- zMMF2
Male2[[3]] <- hMS2
Male2[[4]] <- hMF2
Male2[[5]] <- 0
Nest2 = list()
names(Male2) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel", "Age")
Nest2[[1]] = Female2
Nest2[[2]] = Male2
Nest2[[3]] = SBR
Nest2[[4]] = 0
Nest2[[5]] =  rep(NA, RepAge)
names(Nest2) <- c("Female", "Male", "SBR", "Age", "HelAge")

#F1 <- rep(list(Nest), N)
F1 <- list()
if(SolSocRatio != 1){	
  for(i in c((N*SolSocRatio)+1:(N*(1-SolSocRatio)))){	
    F1[[i]] <- Nest2					
    F1[[i]][[4]] <- sample((RepAge-2):(RepAge+2), 1) #Starting Female	
  } 
}
for(i in c(1:(N*SolSocRatio))){
  print(i)
  F1[[i]] <- Nest1
  F1[[i]][[4]] <- sample((RepAge-2):(RepAge+2), 1) #Starting Female
}
F2 <- list()
x <- Max*DIY*N
for(i in 1:x){
  F2[[i]] <- Nest1
  F2[[i]][[4]] <- sample((RepAge-2):(RepAge+2), 1) #ImmatureFemales
}
F3 <- list()
x <- Max*DIY*N
for(i in 1:x){
  F3[[i]] <- Nest1
  F3[[i]][[4]] <- sample((RepAge-2):(RepAge+2), 1) #MatureFemales
}
M1 <- list()
x <- Max*DIY*N
for(i in 1:x){
  M1[[i]] <- Male1
  M1[[i]][[5]] <- sample((RepAge-2):(RepAge+2), 1) #ImmatureMales
}
M2 <- list()
x <- Max*DIY*N
for(i in 1:x){
  M2[[i]] <- Male1
  M2[[i]][[5]] <- sample((RepAge-2):(RepAge+2), 1) #MatureMales
}
Gen=1
means = c(0, 0, 0)
sds = c(0, 0, 0)
Gens = data.frame(Gen=character(),
                  MeanMSta=character(),
                  MeanMFin=character(), 
                  MeanHSta=character(), 
                  MeanHFin=character(), 
                  SdMSta=character(),
                  SdMFin=character(), 
                  SdHSta=character(), 
                  SdHFin=character(), 
                  stringsAsFactors=FALSE)
MalFem = data.frame(Day=character(),
                    Male=character(),
                    Female=character(), 
                    stringsAsFactors=FALSE)
bs = data.frame(b=character(),
                inc=character(),
                stringsAsFactors = FALSE)


#Run Many Generations
for(i in 1:NumGen){
  if(i == PertNum){
    F1 <- Perturb(F1, MalS02, MalE02, HelS02, HelE02, 0.1)
    print('Yes')
  }  
  Gen = Gen + 1
  if(i%%Skip == 0){
    Data <- Stats(F1)
    Add = data.frame(Gen=Gen, SeasonLength=DIY,MeanMSta = Data[1], MeanMFin = Data[2], MeanHSta = Data[3], MeanHFin = Data[4], SdMSta = Data[5], SdMFin = Data[6], SdHSta = Data[7], SdHFin = Data[8])
    Gens = rbind(Gens, Add) 
  }
  if(Data[3] > 0.05 & Data[3] < 0.95){
    F1 <- NextGeneration(F1)
  }
}
hiveq(Gens, DIY, SolSocRatio, Reps)