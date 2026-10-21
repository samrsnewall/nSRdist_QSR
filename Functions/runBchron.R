#!/usr/bin/env Rscript

## This code is written to be run using the terminal, so that it can be accessed through MATLAB.
## The MATLAB code will invoke the function along with some arguments, which will include:
## - args[1] = the sandbox path
## - args[2] = folder name which has "Inputs" folder and "Outputs" folder
## - args[3] = the core name, so it can access the correct folder in inputs
## - args[4] = the chosen calibration curve
## - args[5] = the desired Delta R (Marine Reservoir Offset)

#The following line is useful for when you want to run this code manually in R, for a given core
#args <- c('/Volumes/ExtDrive850X/MATLAB/nSRdist_code', 'BchronFolders/Bchron_PFandLin_R200M20_Feb4', 'KNR140-51GGC', 'Marine20', '200', '0.1')

## Input corename from command line
args = commandArgs(trailingOnly = TRUE)

##Install the necessary packages
require(Bchron)
require (IntCal)

##Create the marine09 calibration curve and add to the Bchron CalCurve set
marine09 <- IntCal::ccurve(cc = "Marine09",postbomb = FALSE)
createCalCurve("marine09", marine09[,1], marine09[,2], marine09[,3])
file.copy(from = "marine09.rda", to = system.file("data", package = "Bchron"), overwrite = TRUE)

#Define a function to find the maximum value from "densities" list and corresponding age (mode of Bchron output)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Get radiocarbon data
sandboxPath <- args[1L]
folderName <- args[2L]
rdataPath = file.path(sandboxPath, folderName, 'Inputs')
corename <- args[3L]
dataFilename <- paste(corename, '_radiocarbon.txt', sep = "")
rdata = read.delim(file.path(rdataPath, dataFilename))

## Set up directory to store results
genPath <- file.path(sandboxPath, folderName, 'Outputs')
subDir <- corename
dir.create(file.path(genPath, subDir))

## Set up Bchronology Choices
ccurveChosen <- args[4L]
deltaR <- as.numeric(args[5L])
depthSpacing <- as.numeric(args[6L])

##Run Bchronology
chron <- Bchronology(
  ages = rdata$Age, 
  ageSds = sqrt(rdata$Error^2 + deltaR^2),
  positions = rdata$Depth,
  positionThicknesses = rdata$Thickness,
  calCurves = rep(ccurveChosen, length(rdata$Age)),
  outlierProbs = rdata$Outlier1, 
  allowOutside = TRUE, 
  predictPositions = seq(rdata$Depth[1], rdata$Depth[length(rdata$Depth)], by = depthSpacing),
                    )

## Get median ages of calibrated radiocarbon ages

calAges_samples <- sampleAges(chron$calAges)
calAges_medians <- apply(calAges_samples, 2, quantile, prob = c(0.5))
calAges_medians_df <- stack(calAges_medians)[, c("ind", "values")]
names(calAges_medians_df) <- c("Date", "pct50")

##Write useful information out as txt files
write.csv(chron[["thetaPredict"]], file.path(genPath, subDir, "theta.csv"),
          row.names = FALSE)

write.csv(chron[["predictPositions"]], file.path(genPath, subDir, "predictPositions.csv"),
          row.names = FALSE)

write.csv(chron[["phi"]], file.path(genPath, subDir, "phi.csv"),
          row.names = FALSE)

write.table(rdata, file.path(genPath, subDir, "inputData.txt"),
            row.names = FALSE)

write.csv(calAges_medians_df, file.path(genPath, subDir, "calAgesMedian.csv"),
          row.names = FALSE)

