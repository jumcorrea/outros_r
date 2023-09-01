#Sp florestal spp. 
#Progenies de polinizacao aberta ("meios-irmaos")- 5 Anos
# Teste de progenies -
# Delineamento em blocos casualizados (DBC) - 6 blocos
# 27 progenies
# Parcela linear - 5 plantas/parcela
#################################################################
#Pacotes utilizados nas analises 
library(readxl)
library(lmerTest) 
library(dplyr) 
library(plyr)
library(ggplot2)
library(tidyverse)
library(writexl) 

#Para apagar a memoria do R 
rm(list = ls()); ls() 

#Definir diretorio
getwd() 

#Importando os dados 
dados<-read.table("dados2.txt", header = TRUE)
#dados$parcela<-c(dados$Bloco, dados$Prog)

#Estrutura dos dados
dados%>%str()
#Transformando em fatores
dados <- transform(dados,
                   Bloco=as.factor(Bloco),
                   Prog=as.factor(Prog),
                   Ind=as.factor(Ind),
                   Parcela=as.factor(Parcela),
                   Arv=as.factor(Arv), 
                   DAP=as.numeric(DAP))

dados <- na.omit(dados) #deixa s? as plantas vivas

#Estrutura dos dados
dados%>%str()

# Analise grafica dos dados 
  ggplot(dados, aes(Prog,DAP))+
  geom_point(color = "#006633", size = 2, 
             alpha = 0.5)+
  ggtitle("Plot dispersao")+
  labs(x="Progenies",  y="DAP")+
  theme_classic()
  
#Boxplot
ggplot(dados, aes(Prog,DAP))+
  geom_boxplot(alpha = 0.5, 
               outlier.shape  =  NA)+
  ggtitle("Plot dispersao")+
  labs(x="Progenies",  y="DAP")+
  theme_classic()

ggplot(dados, aes(Prog,DAP, fill=Prog))+
  geom_boxplot(alpha = 0.5, 
               outlier.shape  =  NA)+
  labs(x="Progenies",  y="DAP", 
               fill="none")+
  theme(legend.position = "none")

#Verificar os residuos 
ggplot(dados, aes(sample = DAP))+
  geom_qq(color = "#006633", size = 2, 
          alpha = 0.5)+
  geom_qq_line()+
  xlab("Residuos (DAP)")

#Ajuste do modelo misto 
#Prog, Proc e Parcela como efeito aleatorio
m1 <- lmer(DAP ~ Bloco + (1|Prog) +(1|Parcela), dados)
#Resumo do ajuste do modelo 
summary(m1)

#Verificar a significancia dos efeitos aleatorios e fixos
LRT=ranova(m1) 
anova=anova(m1) 

#media fenotipica
media <- mean(dados$DAP, na.rm = T) 

#variance components
vc <- as.data.frame(VarCorr(m1), comp="Variance") 
vc <- subset(vc, select = c("grp","vcov"))
vc

vEntre = vc[1,2] #parcela
vGen = vc[2,2] #genetica
ve= vc[3,2] #erro

#numero de blocos
nb <- length(levels(dados$Bloco)) 
nb

#numero de arvores por parcela
narv <-length(levels(dados$Arv))
narv

va = 4*vGen 
vf = vEntre + va + ve 
h2i = va/vf 
h2d = (0.75*va) / (0.75*va+ve)
h2m = (0.25*va)/(0.25*va + (vEntre/nb) + ve/(narv*nb)) 
acurhi = sqrt(h2i) 
acurhm=sqrt(h2m) 
cvg  = (sqrt(vGen)/media)*100 
cvgi = (sqrt(va)/media)*100 
cve = (sqrt(((0.75*va+ve)/narv)+vEntre))/media*100 
b = cvg/cve # coeficiente de variacao relativo- Vencovsky 1978

###########Parametros geneticos 
Parametros = c("variancia genetica aditiva",
               "Variancia genetica",
               "Variancia residual",
               "Variancia Fenotipica",
               "Herdabilidade no sentido restrito",
               "Herdabilidade dentro de parcela",
               "Herdabilidade media de progenies",
               "Acuracia seletiva (h2i)",
               "Acuracia seletiva (h2m)",
               "Coeficiente de variacao genetico",
               "Coeficiente de variacao genetica aditiva",
               "Coeficiente de variacao experimental",
               "Coeficiente de variacao relativo",
               "Media")           

Valores <- round(c(va,vGen,ve,vf, h2i,h2d,h2m,acurhi,acurhm, cvg,cvgi,cve,b,media), 3) # 3 casas decimais               
parametros_geneticos <- data.frame(Parametros, Valores)
parametros_geneticos

###############################################################
                     # BLUP's - Familias  #  
###############################################################
BLUP = ranef(m1) 
BLUPS = BLUP$Prog
row.names <- data.frame(rownames(BLUPS))
names(row.names) = c("Prog")
BLUP <- data.frame(row.names,g=BLUPS[,1])
BLUP_TRAT <- arrange(BLUP,desc(BLUP$g)) # Ordernar por BLUP
BLUP_TRAT$'u+a'<- BLUP_TRAT$g+media
df_proce= select(dados,Prog)
BLUP = merge(x=df_proce, y=BLUP_TRAT, by="Prog",all.x=T)
BLUP <- BLUP %>% distinct() 
BLUP_TRAT <- arrange(BLUP,desc(BLUP$g)) 


View(BLUP_TRAT) #Blup progenie

############ Plot Ranking Progenies ############
BLUP <- BLUP_TRAT %>%
  top_n(BLUP_TRAT$`u+a`, n=25)

binomnames <- expression(paste("Ranking progenies de", italic(" Especie Abcd"),))

ggplot(BLUP, aes(x= reorder(Prog, -`u+a`), y=`u+a`, fill=Prog))+
  geom_col(alpha=0.5)+
  scale_y_continuous(limits = c(0,21)) +
  labs(fill = "", caption = "DAP - 5 anos")+
  xlab("\nProgenies") +
  ylab("BLUP (DAP)") +
  ggtitle("") +
  labs(title=binomnames)+
  theme(legend.position = "none")


###################### BLUP's Individuais ############################ 
df= select(dados, Bloco, Prog,Arv, DAP)
dfblup = merge(x=df, y=BLUP, by="Prog",all.x=T)
dfblup$'u+a' = dfblup$g+media
dfblup = na.exclude(dfblup)
BLUP_Ind <- arrange(dfblup,desc(dfblup$g)) 
View(BLUP_Ind)

#########################################################
         #GANHOS DE SELECAO 
##########################################################

NIS =100 #Numero de individuos selecionados


#Diferencial de selecao
Ds = mean(BLUP_Ind[1:NIS,6])-media
Ds
#Ganho de selecao individual
Gs = Ds
Gs
#Ganho de selecao percentual
Gs.per = (Gs/media)*100
Gs.per
#Media da populaco melhorada
pm = media + Gs
pm
Ganho <- data.frame(NIS,Gs,Gs.per,pm)
Ganho

##Exportando os resultados 
# output 

write_xlsx(
  list(
    LRT=LRT,
    TesteF=anova,
    parametros_geneticos = parametros_geneticos,
    BLUP_genitores = BLUP_TRAT,
    Ranking_individual = BLUP_Ind),
  "resultados_exemploTP5plantaDAP.xlsx"
)

