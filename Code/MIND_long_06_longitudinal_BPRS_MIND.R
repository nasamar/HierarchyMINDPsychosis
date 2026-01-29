################################################################################
# Script to perform generalized linear modeling (glm) and generalized linear 
# mixed modeling (glmer) on MIND networks and psychiatric symptoms 
################################################################################

#  Copyright (C) 2026 University of Seville
# 
#  Written by Natalia García San Martín (ngarcia1@us.es)
# 
#  This file is part of Hierarchy Longitudinal MIND Gradients Psychosis toolkit.
#
#  Hierarchy Longitudinal MIND Gradients Psychosis toolkit is free software: 
#  you can redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
# 
#  Hierarchy Longitudinal MIND Gradients Psychosis toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Hierarchy Longitudinal MIND Gradients Psychosis toolkit. If not, see 
#  <https://www.gnu.org/licenses/>.

######################  Local structural Longitudinal Analisis ###############################################################


# 1. Initial preparation
####################################################################################

# 1.1. Clean Memory and load libraries
rm(list=ls())

library(lmerTest)        

location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/R/'
source(paste0(location,'regional_brainmap_representation_subcortical.R'))


morphometry <- 'degrees' # do not modify
combat <- 'COMBATLS_covars'

parcellation = 'aparc_500_sym'
# parcellation = 'subcortical'


dx_severity <-'dx'

dx<-'SSD'


x_var <- 'degrees'
x_var <- 'gradients'



gradient <- '_G2'
# gradient <- ''
gradients <- paste0(x_var,gradient)

normalization <- ''
normalization <- '/FEP+CN'
# normalization <- '/FEP_CN'
# normalization <- '/CN'

raw_or_residuals <- 'raw'
# raw_or_residuals <- 'residuals'

sessions <- c(1,2,3,5,10,15,20)
# sessions <- c(1,2,3,10)


cognition_var <- 'BPRS'


GeneralizedLinearModelGauss <- function(D,type){
  if (type=='mixed'){
    names_Objetives <- colnames(D)[c(22:ncol(D))]
    
  }else if(grepl('global',type)){
    names_Objetives <- paste0('global_',x_var)
    
  }else{
    names_Objetives <- colnames(D)[c(22:(ncol(D)-1))] 
  }

  # Names of variables used in interpolation. Used as names of columns in interpolationResults table
  if (!grepl('no mixed',type)){  
    v<-c('Age_inclusion','Sex1','eTIV','mean_euler','CPZ_equivalent','Treatment_Time',x_var,paste0(x_var,':Treatment_Time'),'CPZ_equivalent:Treatment_Time') 
      
    
    if (all(is.na(D[D[,'protocolChange15'],eval(parse(text = cognition_var))]))){
      v<-c(v,'protocol_Change15')
    }
    
  }else{
    v<-c('Age_inclusion','Sex1','eTIV','mean_euler',paste0('global_',x_var))
    
  }
  
  
  # construction of Table for saving beta and p values
  res <- data.frame(matrix(0,ncol=length(names_Objetives),nrow=length(v)))
  names(res)<-names_Objetives
  row.names(res)<-v
  
  resPValues <- res

  n<-length(v)+1
  if (!grepl('global',type)){
    
    for(i in 1:length(names_Objetives)){
      y<-paste0(names_Objetives[[i]])
      D[y]<-scale(D[y])
      if (type=='mixed'){
        centile_or_degree<-paste0('Case:',y)
        formulaText<-paste0(eval(parse(text = cognition_var)),'~1+Age_inclusion+Sex+eTIV+mean_euler+CPZ_equivalent+Treatment_Time+',gsub('Case:','',centile_or_degree),'+',gsub('Case:','',centile_or_degree),':Treatment_Time+CPZ_equivalent:Treatment_Time',ifelse(all(is.na(D[D[,'protocolChange15'],eval(parse(text = cognition_var))])),'+protocolChange15',''),'+(1|Subject)',sep='')
        
        linealModel <- glmer(formulaText,data=D,family=poisson)
          
      }else if (type=='no mixed regional'){
        
        formulaText<-paste0(eval(parse(text = cognition_var)),'~1+Age_inclusion+Sex+eTIV+mean_euler+',y,sep='')
        linealModel <- glm(formulaText,data=D,family=poisson)
          
        
      }
      print(summary(linealModel))
      res[,i]<-as.numeric(coef(summary(linealModel))[2:n,ncol(coef(summary(linealModel)))-1])
      resPValues[,i]<-as.numeric(coef(summary(linealModel))[2:n,ncol(coef(summary(linealModel)))])

    }
  }else{
    D[,paste0('global_',x_var)]<-scale(D[,paste0('global_',x_var)])
    
    if (type=='mixed global'){
      
        formulaText<-paste0(eval(parse(text = cognition_var)),'~1+Age_inclusion+Sex+eTIV+mean_euler+CPZ_equivalent+Treatment_Time+global_',x_var,'+global_',x_var,':Treatment_Time+CPZ_equivalent:Treatment_Time',ifelse(all(is.na(D[D[,'protocolChange15'],eval(parse(text = cognition_var))])),'+protocolChange15',''),'+(1|Subject)',sep='')
        linealModel <- glmer(formulaText,data=D,family = poisson)
      
    }else if (type=='no mixed global'){
      
      formulaText<-paste0(eval(parse(text = cognition_var)),'~1+Age_inclusion+Sex+eTIV+mean_euler+global_',x_var,sep='')
      linealModel <- glm(formulaText,data=D,family=poisson)
      
    }
    print(summary(linealModel))

    res[,1]<-as.numeric(coef(summary(linealModel))[2:n,ncol(coef(summary(linealModel)))-1])
    resPValues[,1]<-as.numeric(coef(summary(linealModel))[2:n,ncol(coef(summary(linealModel)))])
  }
  
  # transpose results tables for correction of p-values
  res<-t(res)
  res_complete <- res
  resPValues<-t(resPValues)

  # FRD correction of P values and selection of significant betas
  for(z in 1:ncol(res)){
    resPValues[,z]<- p.adjust(resPValues[,z], method="fdr")
  }
  
  
  for(i in 1:nrow(res)){
    for(j in 1:ncol(res)){
      if(resPValues[i,j]>0.05){
        res[i,j]<-NaN  # we set to NULL the beta parameters of the variables with p values < 0.05
      }
    }
    
  }
  
  return(list(res = res, res_complete = res_complete))
}


# 1.3. DATA LOADING

Data_Base_covariates <-'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/Covariates_complete.csv'
covariates <-read.csv(Data_Base_covariates,row.names = 1)

Euler <- read.csv('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/euler_longitudinal.csv',row.names = 1)
rownames(Euler) <- gsub('ses-','',gsub('sub-','',rownames(Euler)))

eTIV <- read.csv(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/aparc_formatted_1293s_fs_v7.4.1_long_FUP_13-01-2025/COMBATLS/global_aparc_',combat,'.csv'),row.names = 1)
eTIV <- eTIV[,'EstimatedTotalIntraCranialVol',drop=FALSE]


# Degrees
if (x_var=='degrees'){
  Data_Base_FEP <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/3. MIND_long/Data/degree/',parcellation,'/',combat,'/degree_68_FEP.csv')
  Data_Base_CN <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/3. MIND_long/Data/degree/',parcellation,'/',combat,'/degree_68_CN.csv')
}else if(x_var=='gradients'){  
  Data_Base_FEP <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients_normalized/',parcellation,'/',combat,normalization,'/',raw_or_residuals, '/',gradients,'_FEP.csv')
  Data_Base_CN <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients_normalized/',parcellation,'/',combat,normalization,'/',raw_or_residuals, '/',gradients,'_CN.csv')
}

var_FEP <-read.csv(Data_Base_FEP, row.names = 1)
var_CN <-read.csv(Data_Base_CN, row.names = 1)
morphometry_var <- rbind(var_CN,var_FEP)

covariates <- covariates[rownames(morphometry_var),]
Euler <- Euler[rownames(covariates),]
eTIV <- eTIV[rownames(covariates),]

covariates$Treatment_Time <- covariates$Age_MRI-covariates$Age_inclusion
covariates$eTIV <- eTIV


morphometry_var_longitudinal <- data.frame(morphometry_var[rownames(covariates),],row.names = rownames(covariates))

covariates_longitudinal <- covariates[,c(1:9,which(colnames(covariates)%in%c('eTIV','CPZ_equivalent','Treatment_Time','protocolChange15','Global_Cognitive_Functioning_average','BPRS_Total','SAPS_Total','SANS_Total')))]


Data_Base_Global <- merge(merge(covariates_longitudinal,Euler,by=0),morphometry_var_longitudinal,by.x = 1,by.y = 0)

Data_Base_Global <- Data_Base_Global[Data_Base_Global$Assessment %in% sessions,]
rownames(Data_Base_Global) <- Data_Base_Global[,1]



BPRS <- 'BPRS_Total'

# LONGITUDINAL
{
# for participants in base line we set medication to zero
Data_Base_Global$CPZ_equivalent<-ifelse(Data_Base_Global$Assessment==1,0,Data_Base_Global$CPZ_equivalent)



# To Controls participants we assigned medication  Zero
Data_Base_Global$CPZ_equivalent<-ifelse(Data_Base_Global$Case==0,0.0,Data_Base_Global$CPZ_equivalent)
# fixed dx variable as an R factor
Data_Base_Global$Case<-as.factor(Data_Base_Global$Case)
# Normalize continuous relevant variables
Data_Base_Global$Age_inclusion<-scale(Data_Base_Global$Age_inclusion)
Data_Base_Global$Treatment_Time<-scale(Data_Base_Global$Treatment_Time)
Data_Base_Global$CPZ_equivalent<-scale(Data_Base_Global$CPZ_equivalent)
Data_Base_Global$eTIV<-scale(Data_Base_Global$eTIV)
# fixed Sex variable as an R factor
Data_Base_Global$Sex<-as.factor(Data_Base_Global$Sex)

Data_Base_Global_long <- Data_Base_Global[rownames(covariates_longitudinal),]

if (!(file.exists(paste0(location,"Connectivity/",parcellation,"_longitudinal/",combat,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda")) |
    file.exists(paste0(location,"Connectivity/",parcellation,"_longitudinal/gradients",normalization,'/',raw_or_residuals,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda")))){
  
    InterpolationResults <- GeneralizedLinearModelGauss(Data_Base_Global,'mixed')
    if (x_var == 'degrees'){
      save(InterpolationResults,file= paste0(location,"Connectivity/",parcellation,"_longitudinal/",combat,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda"))
      write.csv(t(InterpolationResults$res_complete),file= paste0(location,"Connectivity/",parcellation,"_longitudinal/",combat,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".csv"))
      write.csv(InterpolationResults,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/degrees/',parcellation,'/InterpolationResults_',gradients,'_',cognition_var,'.csv'),row.names = TRUE)
      
      }else{
      save(InterpolationResults,file= paste0(location,"Connectivity/",parcellation,"_longitudinal/gradients",normalization,'/',raw_or_residuals,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda"))
      write.csv(t(InterpolationResults$res_complete),file= paste0(location,"Connectivity/",parcellation,"_longitudinal/gradients",normalization,'/',raw_or_residuals,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".csv"))
      write.csv(InterpolationResults,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients',ifelse(grep('CN',normalization),'_normalized',''),'/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/InterpolationResults_',gradients,'_',cognition_var,'.csv'),row.names = TRUE)
      
    }
    
    print(InterpolationResults)
    
}else {
  if (x_var == 'degrees'){
    load(paste0(location,"Connectivity/",parcellation,"_longitudinal/",combat,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda"))
  }else{
    load(paste0(location,"Connectivity/",parcellation,"_longitudinal/gradients",normalization,'/',raw_or_residuals,"/InterpolationResults_",x_var,gradient,"_",cognition_var,"_",dx_severity,"_",combat,".Rda"))

  }
  print(InterpolationResults)

}

Data_Base_Global[,paste0('global_',x_var)] = apply(as.matrix(Data_Base_Global[,22:ncol(Data_Base_Global)]), 1, mean, na.rm = TRUE)
InterpolationResults_global <- GeneralizedLinearModelGauss(Data_Base_Global,'mixed global')
print(InterpolationResults_global)


}


# BASELINE
Data_Base_Global_baseline = Data_Base_Global[Data_Base_Global$Assessment==1,]

Data_Base_Global_baseline[,paste0('global_',x_var)] = apply(as.matrix(Data_Base_Global_baseline[,22:ncol(Data_Base_Global_baseline)]), 1, mean, na.rm = TRUE) 

InterpolationResults_baseline_regional <- GeneralizedLinearModelGauss(Data_Base_Global_baseline,'no mixed regional')
if (x_var!='degrees'){
  write.csv(InterpolationResults_baseline_regional,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients',ifelse(grep('CN',normalization),'_normalized',''),'/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/InterpolationResults_',gradients,'_baseline_',cognition_var,'.csv'),row.names = TRUE)
}else{
  write.csv(InterpolationResults_baseline_regional,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/degrees/',parcellation,'/InterpolationResults_',gradients,'_baseline_',cognition_var,'.csv'),row.names = TRUE)
  
}

InterpolationResults_baseline_global <- GeneralizedLinearModelGauss(Data_Base_Global_baseline,'no mixed global')
print(InterpolationResults_baseline_global)



##############################################################################################################
# 3.1 . Representation of data
if (parcellation == 'subcortical'){
  variables<- c('Age_inclusion','Sex1','eTIV',paste0('global_',x_var))
  
  for (i in variables){
    
    data<-data.frame(values=InterpolationResults_baseline_regional$res_complete[,i],sig=!is.na(InterpolationResults_baseline_regional$res[,i]))
    
    if (parcellation=='subcortical'){
      map<-regional_brainmap_representation_subcortical(data,paste0('Baseline ', cognition_var, ' ',dx_severity,' ', x_var,'-> Fixed effect: ',i),max(InterpolationResults_baseline_regional$res_complete,na.rm = TRUE),min(InterpolationResults_baseline_regional$res_complete,na.rm = TRUE),0,'horizontal')
      
    }
    print(map)
    ggsave(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Molecular/plots/Connectivity_longitudinal/Baseline_', cognition_var, ' ',dx_severity,' ' ,x_var,'_Fixed_effect_',gsub(':','_',i),'.png'), width = 8, height = 5, dpi = 300)
    
  }
  
  variables<- c('Age_inclusion','Sex1','eTIV','CPZ_equivalent','Treatment_Time',x_var,paste0(x_var,':Treatment_Time'),'CPZ_equivalent:Treatment_Time',
                paste0(x_var,':CPZ_equivalent'),'protocolChange15')
  
  for (i in variables){
    
    data<-data.frame(values=InterpolationResults$res_complete[,i],sig=!is.na(InterpolationResults$res[,i]))
    
    if (parcellation=='subcortical'){
      map<-regional_brainmap_representation_subcortical(data,paste0('Local ',cognition_var,' ',dx_severity,' vs ',x_var,' -> Fixed effect: ',i),max(InterpolationResults$res_complete),min(InterpolationResults$res_complete),0,'horizontal')
      
    }
    print(map)
    ggsave(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Molecular/plots/Connectivity_longitudinal/Local_', cognition_var,' ',dx_severity,' vs ',x_var,'_Fixed_effect_',gsub(':','_',i),'.png'), width = 8, height = 5, dpi = 300)
    
  }

}
