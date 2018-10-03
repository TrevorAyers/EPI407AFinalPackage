#' ShapeOfExposureAssociation
#'
#' @desciption ShapeOfExposureAssociation is used to find the proper fit of a
#'   continuous exposure variable to a continuous outcome.
#' @param data Data set containing the information being modeled
#' @param Exposure Column name of the dataset being used as the exposure of
#'   interest. This is the variable which will be modified across each of the 5
#'   model types. Example: Exposure = "Exposure1"
#' @param Outcome Column name of the dataset being used as the outcome of
#'   interest. Example: Outcome = "Outcome1"
#' @param ForceModel Options include: "Linear", "Quadratic", "Cubic",
#'   "Log-Linear", and "Quartiles". Sets the Final model equal to the exposure
#'   relationship desired. If ForceModel = NULL the final model will be selected
#'   based on lowest Bayesiant Information Criteria. Default = NULL
#' @param AdjustmentVar Vector of column names used as adjustment variables in
#'   the model. Example: c("adjvar1","adjvar2",...)
#' @param CompareVal Vector of two values of the exposure to compare across the
#'   adjusted models.
#' @details This function creates 5 potential models of the Exposure ~ Outcome
#'   relationship and adjusts for a given set of adjustment variables. The 5
#'   exposure relationships being modeled are: 1. Linear, 2. Quadratic, 3.
#'   Cubic, 4. Log-Linear, and 5. Quartiles of exposure. This function then
#'   selects one of the models and generates a predictive graph of the Exposure
#'   and Outcome relationship and fits a line given the selected model. By
#'   default the appropriate model is chosen based on lowest Bayesian
#'   Information Criteria (BIC).
#' @return
#' \item{BICVals}{List of Bayesian Information Criteria for each of the 5
#' models}
#' \item{Difference}{List of mean change in outcome between the two exposure
#' levels specified in the argument CompareVal for each of the 5 generated
#' models. Only generated if an argument is provided for CompareVal.
#' E(Y|Exp=CompareVal[2])-E(Y|Exp=CompareVal[1]).}
#' \item{FinalModel}{Data Frame containing Estimates, 95 percent confidence
#' intervals, and p-values for each of the covariates of the final model}
#' \item{PlotAssociation}{Scatterplot of Exposure (x) by Outcome (y) with a
#' prediction line fit to the selected model. Generated using ggplot2}
#' @examples
#' ShapeOfExposureAssociation(data=data,Exposure="Age",Outcome="DiastolicBP",
#'     AdjustmentVar=c("Smoked","WalkBike","Gender"),CompareVal=c(20,30))

ShapeOfExposureAssociation<-function(data = NULL, Exposure = NULL, Outcome = NULL,
                                     ForceModel = NULL, AdjustmentVar=NULL, CompareVal = NULL){
  #Create the models
  #Rearrange Adjustment Variables
  if (is.null(AdjustmentVar)){
    adjVars<-NULL
  } else {
    adjVars<-paste(AdjustmentVar,collapse = " + ")
    adjVars<-paste(" + ", adjVars)
  }
  #Linear
  mod1Form<-paste(Outcome, " ~ ", Exposure, adjVars)
  mod1<-lm(as.formula(mod1Form) , data = data)

  #Quadratic
  mod2Form<-paste(Outcome, " ~ ", Exposure, " + ", "I(",Exposure ,"^2)", adjVars)
  mod2<-lm(as.formula(mod2Form) , data = data)

  #Cubic
  mod3Form<-paste(Outcome, " ~ ", Exposure, " + ", "I(",Exposure ,"^2)", " + ", "I(",Exposure ,"^3)", adjVars)
  mod3<-lm(as.formula(mod3Form) , data = data)

  #Log-Linear
  mod4Form<-paste(Outcome, " ~ ", "log(",Exposure,")", adjVars)
  mod4<-lm(as.formula(mod4Form), data = data)

  #Quartiles
  ExposureQuartiles<-quantile(data[,Exposure],probs = c(0.0,0.25,0.5,0.75,1.0),na.rm = TRUE)
  mod5Form<-paste(Outcome, " ~ ", "I(", Exposure,">ExposureQuartiles[2] & ", Exposure, "<=ExposureQuartiles[3]) +
                  I(",Exposure,">ExposureQuartiles[3] & ", Exposure, "<=ExposureQuartiles[4]) +
                  I(",Exposure,">ExposureQuartiles[4] & ", Exposure, "<=ExposureQuartiles[5])",
                  adjVars)
  mod5<-lm(as.formula(mod5Form), data = data)

  #Compare the models
  BICVec<-round(c(BIC(mod1),BIC(mod2),BIC(mod3),BIC(mod4),BIC(mod5)),2)
  tempOut<-data.frame(BICVec)
  colnames(tempOut)<-c("BIC")
  rownames(tempOut)<-c("Linear","Quadratic","Cubic","Log-Linear","Quartiles")

  #Compare estimates
  if (is.null(CompareVal)==FALSE) {
    linAdj  <-(summary(mod1)$coef[2,1]*CompareVal[2])-
      (summary(mod1)$coef[2,1]*CompareVal[1])
    sqrdAdj <-(summary(mod2)$coef[2,1]*CompareVal[2]+summary(mod2)$coef[3,1]*CompareVal[2]^2)-
      (summary(mod2)$coef[2,1]*CompareVal[1]+summary(mod2)$coef[3,1]*CompareVal[1]^2)
    cubeAdj <-(summary(mod3)$coef[2,1]*CompareVal[2]+summary(mod3)$coef[3,1]*CompareVal[2]^2+summary(mod3)$coef[4,1]*CompareVal[2]^3)-
      (summary(mod3)$coef[2,1]*CompareVal[1]+summary(mod3)$coef[3,1]*CompareVal[1]^2+summary(mod3)$coef[4,1]*CompareVal[1]^3)
    logAdj  <-(summary(mod4)$coef[2,1]*log(CompareVal[2]))-
      (summary(mod4)$coef[2,1]*log(CompareVal[1]))
    quartAdj<-(summary(mod5)$coef[2,1]*I(CompareVal[2]>ExposureQuartiles[2] & CompareVal[2]<=ExposureQuartiles[3])+
                 summary(mod5)$coef[3,1]*I(CompareVal[2]>ExposureQuartiles[3] & CompareVal[2]<=ExposureQuartiles[4])+
                 summary(mod5)$coef[4,1]*I(CompareVal[2]>ExposureQuartiles[4] & CompareVal[2]<=ExposureQuartiles[5]))-
      (summary(mod5)$coef[2,1]*I(CompareVal[1]>ExposureQuartiles[2] & CompareVal[1]<=ExposureQuartiles[3])+
         summary(mod5)$coef[3,1]*I(CompareVal[1]>ExposureQuartiles[3] & CompareVal[1]<=ExposureQuartiles[4])+
         summary(mod5)$coef[4,1]*I(CompareVal[1]>ExposureQuartiles[4] & CompareVal[1]<=ExposureQuartiles[5]))
    tempDiff<-data.frame(c(linAdj,sqrdAdj,cubeAdj,logAdj,quartAdj),row.names = c("Linear","Quadratic","Cubic","Log-Linear","Quartiles"))
    colnames(tempDiff)<-c("AdjustedDiff")
  }

  #Select the correct model
  #Set "Result"
  if (is.null(ForceModel)){
    result<-which.min(BICVec)
  } else if (ForceModel=="Linear"){
    result<-1
  } else if (ForceModel=="Quadratic"){
    result<-2
  } else if (ForceModel=="Cubic"){
    result<-3
  } else if (ForceModel=="Log-Linear"){
    result<-4
  } else if (ForceModel=="Quartiles"){
    result<-5
  }
  #Identify correct model using "Result"
  if (result==1){
    FinalModel<-mod1
  } else if (result==2){
    FinalModel<-mod2
  } else if (result==3){
    FinalModel<-mod3
  } else if (result==4){
    FinalModel<-mod4
  } else if (result==5){
    FinalModel<-mod5
  }
  #Create paste output of Final Model
  finalSum<-summary(FinalModel)$coef
  finalConf<-NULL
  for (i in 1:nrow(finalSum)){
    finalConf[i]<-paste("(",round(confint(FinalModel)[i,1],2),",",round(confint(FinalModel)[i,2],2),")")
  }
  finalOut<-data.frame(cbind(round(finalSum[,1],2),finalConf,finalSum[,4]))
  colnames(finalOut)<-c("Estimate","Estimate95CI","Pval")

  #Plot Association
  seqExposure<-seq(from=min(data[,Exposure],na.rm=TRUE),
                   to=max(data[,Exposure],na.rm=TRUE) ,length.out =1000)
  inputframe<-data.frame(Exposure=seqExposure)

  if (is.null(AdjustmentVar)==FALSE){
    for (i in 1:length(AdjustmentVar)){
      if (is.numeric(data[,AdjustmentVar[i]])==FALSE){
        data[,AdjustmentVar[i]]<-as.character(data[,AdjustmentVar[i]])
        mostfreq<-names(sort(table(data[,AdjustmentVar[i]]),decreasing=TRUE))[1]
        inputframe<-cbind(inputframe,mostfreq)
      } else {
        mediancont<-median(data[,AdjustmentVar[i]],na.rm=TRUE)
        inputframe<-cbind(inputframe,mediancont)
      }
    }
    colnames(inputframe)<-c(Exposure,AdjustmentVar)
  } else {
    colnames(inputframe)<-c(Exposure)
  }

  predictedOutcome<-predict(FinalModel,inputframe)
  predictFrame<-data.frame(Exp=seqExposure,Out=predictedOutcome)

  if (result==5){
    Quart1<-predictFrame[which(predictFrame[,1]>=ExposureQuartiles[1] &
                                 predictFrame[,1]<=ExposureQuartiles[2]),]
    Quart2<-predictFrame[which(predictFrame[,1]>ExposureQuartiles[2] &
                                 predictFrame[,1]<=ExposureQuartiles[3]),]
    Quart3<-predictFrame[which(predictFrame[,1]>ExposureQuartiles[3] &
                                 predictFrame[,1]<=ExposureQuartiles[4]),]
    Quart4<-predictFrame[which(predictFrame[,1]>ExposureQuartiles[4] &
                                 predictFrame[,1]<=ExposureQuartiles[5]),]

    finalPlot<-ggplot(data, aes(data[,Exposure],data[,Outcome]))+
      geom_point(alpha=0.8,col="light blue")+
      geom_line(data = Quart1, aes(Quart1[,1],Quart1[,2]),col="blue",size=2)+
      geom_line(data = Quart2, aes(Quart2[,1],Quart2[,2]),col="blue",size=2)+
      geom_line(data = Quart3, aes(Quart3[,1],Quart3[,2]),col="blue",size=2)+
      geom_line(data = Quart4, aes(Quart4[,1],Quart4[,2]),col="blue",size=2)+
      xlab("Exposure")+
      ylab("Outcome")+
      theme_bw()+
      theme(axis.title = element_text(face = "bold"))
  } else {
    finalPlot<-ggplot(data, aes(data[,Exposure],data[,Outcome]))+
      geom_point(alpha=0.8,col="light blue")+
      geom_line(data = predictFrame, aes(predictFrame[,1],predictFrame[,2]),col="blue",size=2)+
      xlab("Exposure")+
      ylab("Outcome")+
      theme_bw()+
      theme(axis.title = element_text(face = "bold"))
  }
  list(BICVals=tempOut, Difference=tempDiff, FinalModel=finalOut, PlotAssociation=finalPlot)
}

#' Standard Attack Roll
#'
#' @description I created a function which performs the basic decision tree used
#'   when making an attack roll in Dungeons and Dragons 5th Edition.
#'
#' @param Character a data frame containing character information for a Dungeons
#'   & Dragons character. This function specifically referrs to columns "Str"
#'   and "Dex" of this character sheet data frame.
#' @param Range Categorical input for the type of attack being made ("Melee" or
#'   "Ranged")
#' @param profBonus Integer value of the proficiency bonus of the character
#'   making the attack
#' @param Proficiency "Yes" or "No" Is the character proficient with the weapon
#'   they are making the attack with? Default: "Yes"
#' @param Circumstance Refers to circumstances in which the roll is being made.
#'   If the player has the advantage on the roll, Cirumstance="Advantage". If
#'   the player has a disadvantage on the roll, Cirumstance="Disdvantage". For a
#'   standard attack without either, Circumstance=NULL. Default=NULL
#' @param Finesse "Yes" or "No" If the weapon is melee, does it have the finesse
#'   property?
#' @details This function simulates the process of rolling a 20-sided die and
#'   adds the proper stat (Strength or Dexterity depending on weapon the attack
#'   was made with) and proficiency bonuses relevant when making an attack roll
#'   as described in the Dungeons and Dragons 5th Edition Players Handbook (pg
#'   194). This function does not take into account penalties associated with
#'   target range, in its current version.
#' @return A string of text indicating the outcome of the roll before and after
#'   bonuses were added.
#' @examples
#' sampleCharacter<-data.frame(cbind(Name="Trevor",Class="Barbarian",Str=16,Dex=13,Con=14,Int=8,Wis=9,Cha=11))
#' standardAttack(Character = sampleCharacter, Range = "Melee", Circumstance =
#' "Advantage", Finesse = "No", profBonus = 2)

standardAttack<-function(Character=NULL,Range=NULL,profBonus=NULL,Circumstance=NULL,Finesse=NULL,Proficiency="Yes"){
  if (is.null(Circumstance)){
    roll<-ceiling(runif(1, min = 0, max = 20))
  } else if (Circumstance=="Disadvantage") {
    roll<-min(ceiling(runif(2, min = 0, max = 20)))
  } else if (Circumstance=="Advantage"){
    roll<-max(ceiling(runif(2, min = 0, max = 20)))
  }
  if (roll==1){
    print("Critical Failure! See GM for consequences...")
  } else if (roll==20){
    print("Oh baby a critical! Let the crit dice roll!")
  } else {
    print(paste("I rolled a natural...",roll,", which after adding my modifiers becomes..."))
    if (Range=="Melee" & Finesse=="No"){
      statMod<-floor((as.numeric(as.character(Character$Str))-10)/2)
    } else if (Range=="Melee" & Finesse=="Yes" | Range=="Ranged"){
      statMod<-floor((as.numeric(as.character(Character$Dex))-10)/2)
    }
    ifelse(Proficiency=="Yes",
           roll<-roll+profBonus,
           roll<-roll)
    roll<-roll+statMod
    print(roll)
  }
}
