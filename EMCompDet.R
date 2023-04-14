library(plyr)
library(reshape2)
library(ggnewscale)
library(tidyverse)

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

NextGeneration <- function(StartingReps){
  Son <- Male1
  Daughter <- Female1
  NImmRep <- 1
  NMatRep <- as.numeric(length(StartingReps)) + 1
  NImmMal <- 1
  NMatMal <- 1
  NH <- 1
  df = data.frame(matrix(ncol = 4, nrow = DIY))
  x <- c("Day", "Male", "Help", "Repr")
  colnames(df) <- x
  df$Day <- c(1:DIY)
  MatRep[1:(NMatRep-1)] <- StartingReps 
  for(j in 1:DIY){
    for(i in 1:(NMatRep-1)){
      Kids = rpois(1, MatRep[[i]][["SBR"]])
      if(Kids > 0){
        for(p in 1:Kids){
          z1 <- 0.5 * (MatRep[[i]][["Female"]][["Starting Male"]][1] + MatRep[[i]][["Female"]][["Starting Male"]][2]) #Average of male (chrom 1) with male (chrom 2) starting ratios
          z2 <- 0.5 * (MatRep[[i]][["Female"]][["Finishing Male"]][1] + MatRep[[i]][["Female"]][["Finishing Male"]][2]) #Average of male (chrom 1) with male (chrom 2) finishing ratios
          z <- (as.numeric(j))*(z2) + z1
          if(j==72){
            print(z)
          }
          if(z<0){ #Can't have values less than 0 or greater than 1
            z <- 0
          }
          if(z>1){
            z <- 1
          }
          if(runif(1) < z){ #If Male (remember z is male rate)
            if(LD == 0){
              SonMSta <- MatRep[[i]][["Female"]][["Starting Male"]][sample(1:2, 1)] #Pick a Random Chromosome from Female (Haplodiploid!)
              SonMFin <- MatRep[[i]][["Female"]][["Finishing Male"]][sample(1:2, 1)] #Same with the Finishing Prop Rate
              SonHSta <- MatRep[[i]][["Female"]][["Starting Hel"]][sample(1:2, 1)] #Same with the Helping Rate
              SonHFin <- MatRep[[i]][["Female"]][["Finishing Hel"]][sample(1:2, 1)] #Same with the Helping Rate
            }
            else{
              Chr <- sample(1:2, 1)
              SonMSta <- MatRep[[i]][["Female"]][["Starting Male"]][Chr] #Pick a Random Chromosome from Female (Haplodiploid!)
              SonMFin <- MatRep[[i]][["Female"]][["Finishing Male"]][Chr] #Same with the Finishing Prop Rate
              SonHSta <- MatRep[[i]][["Female"]][["Starting Hel"]][Chr] #Same with the Helping Rate
              SonHFin <- MatRep[[i]][["Female"]][["Finishing Hel"]][Chr] #Same with the Helping Rate
            }
            Son[["Starting Male"]] <- mutate(SonMSta) #MutateTheValuesAndCreateASonWithThem
            Son[["Finishing Male"]] <- mutate2(SonMFin)
            Son[["Starting Hel"]] <- mutate(SonHSta)
            Son[["Finishing Hel"]] <- mutate2(SonHFin)
            Son[["Age"]] <- 0
            M1[[NImmMal]] <- Son
            NImmMal <- NImmMal + 1
            #cat("A Male was Born on Day", j)
          }
          else{
            if(LD == 0){
              DaughterMSta <- MatRep[[i]][["Female"]][["Starting Male"]][sample(1:2, 1)] #Pick a Random Chromosome from Female (Haplodiploid!)
              DaughterMFin <- MatRep[[i]][["Female"]][["Finishing Male"]][sample(1:2, 1)] #Same with the Finishing Prop Rate
              DaughterHSta <- MatRep[[i]][["Female"]][["Starting Hel"]][sample(1:2, 1)] #Same with the Helping Rate
              DaughterHFin <- MatRep[[i]][["Female"]][["Finishing Hel"]][sample(1:2, 1)] #Same with the Helping Rate
            }
            else{
              Chr <- sample(1:2, 1)
              DaughterMSta <- MatRep[[i]][["Female"]][["Starting Male"]][Chr] #Pick a Random Chromosome from Female (Haplodiploid!)
              DaughterMFin <- MatRep[[i]][["Female"]][["Finishing Male"]][Chr] #Same with the Finishing Prop Rate
              DaughterHSta <- MatRep[[i]][["Female"]][["Starting Hel"]][Chr] #Same with the Helping Rate
              DaughterHFin <- MatRep[[i]][["Female"]][["Finishing Hel"]][Chr] #Same with the Helping Rate
            }
            Daughter[["Starting Male"]][1] <- mutate(DaughterMSta) #MutateTheValuesAndCreateASonWithThem
            Daughter[["Finishing Male"]][1] <- mutate2(DaughterMFin)
            Daughter[["Starting Hel"]][1] <- mutate(DaughterHSta)
            Daughter[["Finishing Hel"]][1] <- mutate2(DaughterHFin)
            DaughterMSta <- MatRep[[i]][["Male"]][["Starting Male"]]
            DaughterMFin <- MatRep[[i]][["Male"]][["Finishing Male"]]
            DaughterHSta <- MatRep[[i]][["Male"]][["Starting Hel"]]
            DaughterHFin <- MatRep[[i]][["Male"]][["Finishing Hel"]]
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
                if (is.na(MatRep[[i]][["HelAge"]][k])){
                  MatRep[[i]][["HelAge"]][k] = 0
                  #cat("A Helper was Born on Day", j)
                  break
                }
              }
            }
            else{
              ImmRep[[NImmRep]][["Female"]] <- Daughter
              ImmRep[[NImmRep]][["SBR"]] <- BirthRate
              ImmRep[[NImmRep]][["HelAge"]] <- rep(NA, RepAge)
              ImmRep[[NImmRep]][["Age"]] <- 0
              NImmRep = NImmRep + 1
              #cat("A Reproductive was Born on Day", j)
            }
          }
        }
      }
      if(runif(1) < Sfe){
        MatRep[[i]] = MatRep[[NMatRep-1]]
        NMatRep = NMatRep - 1
      }
      else{
        MatRep[[i]][["Age"]] <- MatRep[[i]][["Age"]] + 1
        MatRep[[i]][["HelAge"]] = MatRep[[i]][["HelAge"]] + 1/RepAge
        X = which(MatRep[[i]][["HelAge"]] >= 1) 
        MatRep[[i]][["HelAge"]][X] = NA
        MatRep[[i]][["SBR"]] <- MatRep[[i]][["SBR"]] + length(X)*b
        NH = NH + length(X)
      }
    }
    if(NImmRep > 1){
      for(i in 1:(NImmRep-1)){
        if(ImmRep[[i]][["Age"]] == RepAge + TimeToNest){
          if(1==1){
            #cat("A Reproductive Matured on Day", j)
            MatRep[[NMatRep]][["Female"]] <- ImmRep[[i]][["Female"]]
            MatRep[[NMatRep]][["Male"]] <- M2[[sample(1:NMatMal, 1)]]
            NMatRep = NMatRep + 1
            ImmRep[[i]] = ImmRep[[NImmRep-1]]
            ImmRep[[i]][["Age"]] = ImmRep[[i]][["Age"]] + 1
            NImmRep = NImmRep - 1
          }
          else{
            ImmRep[[i]][["Age"]] = ImmRep[[i]][["Age"]] - 2 #If no males are mature, hold off
          }
        }   #Fertilize Mature Reproductive Females
        else{
          if(runif(1) < Sfe){
            ImmRep[[i]] = ImmRep[[NImmRep-1]]
            NImmRep = NImmRep - 1
            ImmRep[[i]][["Age"]] = ImmRep[[i]][["Age"]] + 1
          }
          else{
            ImmRep[[i]][["Age"]] = ImmRep[[i]][["Age"]] + 1
          }
        }
      }
    }
    
    if(NImmMal > 1){
      for(i in 1:(NImmMal-1)){
        if (M1[[i]][["Age"]] == RepAge){
          M2[[NMatMal]] = M1[[i]] 
          NMatMal = NMatMal + 1
          M1[[i]] = M1[[NImmMal-1]]
          NImmMal = NImmMal - 1
          #cat("A Male Matured on Day", j)
        }
        if(runif(1) < Sme){
          M1[[i]] = M1[[NImmMal-1]]
          NImmMal = NImmMal - 1
          M1[[i]][["Age"]] = M1[[i]][["Age"]] + 1
        }
        else{
          M1[[i]][["Age"]] = M1[[i]][["Age"]] + 1
        }
      }
    }
    if(NMatMal > 1){
      for(i in 1:(NMatMal-1)){
        if(runif(1) < Sme){
          M2[[i]] = M2[[NMatMal-1]]
          NMatMal = NMatMal - 1
        }
      }
    }
    df$Male[j] <- NMatMal-1
    df$Help[j] <- NH-1
    df$Repr[j] <- NMatRep-1
    #df$ImmM[j] <- NImmMal-1
    #df$ImmF[j] <- NImmRep-1
  }
  for(i in 1:N){
    RF =  as.numeric(sample(1:(NMatRep-1), 1))
    StartingReps[[i]] = MatRep[[RF]]
    StartingReps[[i]][["SBR"]] <- BirthRate
    StartingReps[[i]][["HelAge"]] <- rep(NA, RepAge)
    NMatRep = NMatRep-1
    MatRep[[RF]] <- MatRep[[NMatRep]]
  }
  return(df)
}

#This Function will convert set of nests StartingReps to a new generation

Perturb <- function(StartingReps, Perturbation, Percentage){
  Intruders <- floor(length(StartingReps)*Percentage)
  for(i in 1:Intruders){
    StartingReps[[i]][["Female"]][["Starting Hel"]][1] <- StartingReps[[i]][["Female"]][["Starting Hel"]][1] + Perturbation
    StartingReps[[i]][["Female"]][["Starting Hel"]][2] <- StartingReps[[i]][["Female"]][["Starting Hel"]][2] + Perturbation
    StartingReps[[i]][["Male"]][["Starting Hel"]] <- StartingReps[[i]][["Male"]][["Starting Hel"]] + Perturbation
  }
  return(StartingReps)
}

#Set Variables
b = 0.15 #Increase In Summer Growth Rate from Helpers
BirthRate = 0.05
Sf = 0 #Proportion of Foundresses that Die before end of summer
Sm = 0 #Mort Rate of Males
DIY <- 120  #Length from diapause > diapause
Sfe = Sf / 120 #Per Day Mort Rate
Sme = Sm / 120 #Per Day Mort Rate
mu = 0 #Mut Rate Per Allele Per Generation
sdmu = 0 #St. Dev. of Mutation Size
sdmu2 = 0
Max = fib(ceiling((DIY) / (1/BirthRate)) +1) #MaximumNumberInABrood
N = 50; #Nests
RepAge <- 27
TimeToNest <- 10
NumGen = 1 #Number Gens
Lag = 0 #Time before helper starts evolving
Skip = 1 #Output Every 10th
MalS01 = 0.712
MalE01 = -0.00548
HelS01 = 0.318
HelE01 = -0.0078
#MalS01 = 0.25
#MalE01 = 0.0025
#HelS01 = 0.66
#HelE01 = -0.008
#MalS01 = 0.0537
#MalE01 = 0.00307
#HelS01 = 0.954
#HelE01 = -0.00750
zFMS1 = c(MalS01, MalS01)
zFMStartingReps = c(MalE01, MalE01)
zMMS1 = MalS01
zMMStartingReps = MalE01
hFS1 = c(HelS01, HelS01)
hFStartingReps = c(HelE01, HelE01)
hMS1 = HelS01
hMStartingReps = HelE01
SolSocRatio = 1
SBR = BirthRate
PertNum = 1000000
LD = 0
ChangeInClimate = 2332
DaysInYearSd = 19991

#Create Data Frame
Female1 = list()
Female1[[1]] <- zFMS1
Female1[[2]] <- zFMStartingReps
Female1[[3]] <- hFS1
Female1[[4]] <- hFStartingReps
names(Female1) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel")
Male1 = list()
Male1[[1]] <- zMMS1
Male1[[2]] <- zMMStartingReps
Male1[[3]] <- hMS1
Male1[[4]] <- hMStartingReps
Male1[[5]] <- 0
Nest1 = list()
names(Male1) <- c("Starting Male", "Finishing Male", "Starting Hel", "Finishing Hel", "Age")
Nest1[[1]] = Female1
Nest1[[2]] = Male1
Nest1[[3]] = SBR
Nest1[[4]] = 0
Nest1[[5]] =  rep(NA, RepAge)
names(Nest1) <- c("Female", "Male", "SBR", "Age", "HelAge")
N1<-N*SolSocRatio
StartingReps <- list()
for(i in c(1:N1)){
  StartingReps[[i]] <- Nest1
  StartingReps[[i]][[4]] <- sample((RepAge*0):(RepAge*0), 1) #Starting Female
}
ImmRep <- list()
x <- Max*DIY*N
for(i in 1:x){
  ImmRep[[i]] <- Nest1
  ImmRep[[i]][[4]] <- sample((RepAge*0):(RepAge*0), 1) #ImmatureFemales
}
MatRep <- list()
x <- Max*DIY*N
for(i in 1:x){
  MatRep[[i]] <- Nest1
  MatRep[[i]][[4]] <- sample((RepAge*0):(RepAge*0), 1)#MatureFemales
}
M1 <- list()
x <- Max*DIY*N
for(i in 1:x){
  M1[[i]] <- Male1
  M1[[i]][[5]] <- sample((RepAge*0):(RepAge*0), 1) #ImmatureMales
}
M2 <- list()
x <- Max*DIY*N
for(i in 1:x){
  M2[[i]] <- Male1
  M2[[i]][[5]] <- sample((RepAge*0):(RepAge*0), 1) #MatureMales
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
for(j in 1:1){
  StartingReps <- list()
  Gen=0
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
  for(i in c(1:N1)){
    StartingReps[[i]] <- Nest1
    StartingReps[[i]][[4]] <- sample((RepAge*0):(RepAge*0), 1) #Starting Female
  }
}
df <- list(NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps), NextGeneration(StartingReps))
means = data.frame(aaply(laply(df, as.matrix), c(2, 3), mean))
means$Male = means$Male / 50
means$Help = means$Help / 50
means$Repr = means$Repr / 50
sd = data.frame(aaply(laply(df, as.matrix), c(2, 3), var))
sd$Male = sqrt(sd$Male) / 50
sd$Help = sqrt(sd$Help) / 50
sd$Repr = sqrt(sd$Repr) / 50
mya <- melt(means, id=c("Day"))
mya$Behavior = "Social"
mya$Model = "Autonomous"
mysd <- melt(sd, id=c("Day"))
mya$sd = mysd$value


#Deterministic
TauM = 28
TauH = 28
TauR = 38
GammaM = 0
GammaR = 0
Alpha = BirthRate
Beta = b
End = DIY
GammaMe = GammaM / 120
GammaRe = GammaR / 120
Start = 50
M0 = 0.712
MPrime = -0.00548
H0 = 0.318
HPrime = -0.0078
#M0 = 0.25
#MPrime = 0.0025
#H0 = 0.66
#HPrime = -0.008
#M0 = 0.0537
#MPrime = 0.00307
#H0 = 0.954
#HPrime = -0.00750
df = data.frame(matrix(ncol = 4, nrow = End+1))
x <- c("Day", "Male", "Help", "Repr")
colnames(df) <- x
SolSocSim <- df
df$Day <- c(0:End)
df$Male[1] = 0
df$Help[1] = 0
#df$ImmM[1] = 0
#df$ImmF[1] = 0
df$Repr[1] = Start
for(i in 1:End){
  if(i>=TauM){
    Pm = M0+(MPrime*(i+1-TauM))
    Pm = min(max(Pm,0),1)
    print(Pm)
    Ph = (1-Pm)*(H0+(HPrime*(i+1-TauH)))
    Ph = min(max(Ph,0),1)
    Pr = 1-Pm-Ph
    df$Male[i+1] = (1-GammaMe)*(df$Male[i]) + ((1-GammaMe)^TauM)*((Alpha*df$Repr[i+1-TauM] + Beta*(df$Help[i+1-TauM]))*Pm)
    df$Help[i+1] = (1-GammaRe)*(df$Help[i]) + ((1-GammaRe)^TauH)*((Alpha*df$Repr[i+1-TauH] + Beta*(df$Help[i+1-TauH]))*Ph)
    #df$ImmM[i+1] = (1-GammaMe)*(df$ImmM[i]) + (Alpha*df$Repr[i] + Beta*(df$Help[i]))*Pm - ((1-GammaMe)^TauM)*((Alpha*df$Repr[i+1-TauM] + Beta*(df$Help[i+1-TauM]))*Pm)
  }
  else{
    Pm = M0+(MPrime*(i))
    Pm = min(max(Pm,0),1)
    Ph = (1-Pm)*(H0+(HPrime*(i)))
    Ph = min(max(Ph,0),1)
    Pr = 1-Pm-Ph
    df$Male[i+1] = (1-GammaMe)*df$Male[i]
    df$Help[i+1] = (1-GammaRe)*df$Help[i]
    #df$ImmM[i+1] = (1-GammaMe)*(df$ImmM[i]) + (Alpha*df$Repr[i] + Beta*(df$Help[i]))*Pm
    
  }
  if(i>=TauR){
    Pm = M0+(MPrime*(i+1-TauR))
    Pm = min(max(Pm,0),1)
    Ph = (1-Pm)*(H0+(HPrime*(i+1-TauR)))
    Ph = min(max(Ph,0),1)
    Pr = 1-Pm-Ph
    df$Repr[i+1] = ((1-GammaRe)*df$Repr[i]) + (((1-GammaRe)^TauR)*((Alpha*df$Repr[i+1-TauR] + Beta*(df$Help[i+1-TauR])))*Pr)
    #df$ImmF[i+1] = (1-GammaRe)*(df$ImmF[i]) + (Alpha*df$Repr[i] + Beta*(df$Help[i]))*Pr - ((1-GammaRe)^TauR)*((Alpha*df$Repr[i+1-TauR] + Beta*(df$Help[i+1-TauR]))*Pr)
  }
  else{
    df$Repr[i+1] = (1-GammaRe)*df$Repr[i]
    #df$ImmF[i+1] = (1-GammaRe)*(df$ImmF[i]) + (Alpha*df$Repr[i] + Beta*(df$Help[i]))*Pr
  }
}
df$Help <- df$Help / 50
df$Male <- df$Male / 50
df$Repr <- df$Repr / 50
myd <- melt(df, id="Day")
myd$Behavior <- "Social"
myd$Model <- "Deterministic"
myd$sd = 0

AutDetData <- rbind(myd, mya)
cols <- c("Male" = "#88cde9", "Help" = "#74e841", "Repr"="#a022b0", "ImmM" = "grey", "ImmF" = "maroon")
p <- ggplot(myd, aes(x=Day, y=value, col=variable, linetype=Model)) + 
  geom_line(aes(col=variable, size='1')) + 
  geom_linerange(data = mya, aes(ymin=value-sd,  ymax=value+sd, colour=variable), alpha=0.3) +
  scale_colour_manual(values = cols) 
p + xlab("Day") + ylab("Number of Offspring") + xlim(0,180) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none', axis.text=element_text(size=14),axis.title=element_text(size=20))
