#Sp florestal spp. 
#Progenies de polinizacao aberta ("meios-irmaos")- 5 Anos
# Teste de  progenies
# Delineamento em blocos casualizados (DBC) - 6 blocos
# 27 progenies 
# 1 planta/parcela
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
dados<-read.table("dados1.txt", header = TRUE)
head(dados)
View(dados)
tail(dados)

#Estrutura dos dados
dados%>%str() 

#Transformando em fatores
dados$Ind=as.factor(dados$Ind)
dados$Bloco=as.factor(dados$Bloco)
dados$Prog=as.factor(dados$Prog)
dados$Arv=as.factor(dados$Arv)
dados$DAP=as.numeric(dados$DAP)

#Estrutura dos dados
dados%>%str()

# Analise grafica dos dados 
ggplot(dados, aes(Prog,DAP, color= Pop))+
  geom_point(color = "#006633", 
             size = 2, alpha = 0.5)+
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
  labs(x="Progenies",  y="DAP", fill="none")+
  theme(legend.position = "none")

#Verificar os residuos 
ggplot(dados, aes(sample = DAP))+
  geom_qq(color = "#006633", size = 2, 
          alpha = 0.5)+
  geom_qq_line()+
  xlab("Residuos (DAP)")

#Normalidade 
shapiro.test(dados$DAP)

#Ajuste do modelo misto 
#Prog, Proc e Parcela como efeito aleatorio
#m1<- lmer(variavel ~ fator fixo + (1|fator aleatorio), conjunto de dados)
m1 <- lmer(DAP ~ Bloco + (1|Prog), dados) 

#Resumo do ajuste do modelo 
summary(m1)

#Verificar a significancia dos efeitos aleatorios e fixos
LRT=ranova(m1) #LRT 
anova=anova(m1) #Teste F 

#media fenotipica
media <- mean(dados$DAP, na.rm = T) 

#variance components
vc <- as.data.frame(VarCorr(m1), comp="Variance") 
vc <- subset(vc, select = c("grp","vcov"))
vc

#vGen = variancia da progenie; ve = variancia do erro (residual)
vGen = vc[1,2] #linha 1, coluna 2
ve= vc[2,2] #linha 2, coluna 2

#numero de blocos
nb <- length(levels(dados$Bloco)) 
nb

va = 4*vGen #variancia genetica (aditiva)
vf = va + ve #variancia fenotipica
h2i = va/vf #h2 individual sentido restrito
h2d = (0.75*va) / (0.75*va+ve) #h2 dentro
h2m = (0.25*va)/(0.25*va + ve/(nb))  #h2 media de progenies
acurhi = sqrt(h2i) #acuracia individual progenies
acurhm=sqrt(h2m)  #acuracia media progenies
cvg  = (sqrt(vGen)/media)*100 #coef variacao genetica
cvgi = (sqrt(va)/media)*100  #coef variacao genetica aditivo
cve = (sqrt(ve)/media)*100 #coef variacao erro
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
BLUP = ranef(m1) # extrair os BLUPS
BLUPS = BLUP$Prog
row.names <- data.frame(rownames(BLUPS))
names(row.names) = c("Prog")
BLUP <- data.frame(row.names,g=BLUPS[,1])                                  
BLUP_TRAT <- arrange(BLUP,desc(BLUP$g)) 
BLUP_TRAT$'u+a'<- BLUP_TRAT$g+media
df_proce= select(dados,Prog)


BLUP = merge(x=df_proce, y=BLUP_TRAT, by="Prog",all.x=T)
BLUP <- BLUP %>% distinct() 
BLUP_TRAT <- arrange(BLUP,desc(BLUP$g))

View(BLUP_TRAT) #Blup progenie

#---- plot ranking progenies ---------
BLUP <- BLUP_TRAT %>%
  top_n(BLUP_TRAT$`u+a`, n=25)

binomnames <- expression(paste("Ranking progenies de", italic(" Especie Abcd"),))

ggplot(BLUP, aes(x= reorder(Prog, -`u+a`), y=`u+a`))+
  geom_col(aes(fill=Prog), alpha=0.5)+
  scale_y_continuous(limits = c(0,60)) +
  labs(fill = "", caption = "DAP - 5 anos")+
  xlab("\nProgenies") +
  ylab("BLUP (DAP)") +
  ggtitle("") +
  labs(title=binomnames)+
  theme(legend.position = "none") 

###################### BLUP's Individuais ############################ 
df= select(dados, Bloco, Prog,DAP)
dfblup = merge(x=df, y=BLUP, by="Prog",all.x=T)
dfblup$yijblup = dfblup$DAP-dfblup$g
dfblup$yijblup = dfblup$yijblup*h2d
dfblup$yijblup = dfblup$yijblup+dfblup$g
dfblup$yijblup = dfblup$yijblup - mean(dfblup$yijblup, na.rm = T)
dfblup = na.exclude(dfblup)
BLUP_Ind <- arrange(dfblup,desc(dfblup$yijblup)) 

BLUP_Ind <- BLUP_Ind[,-6]

View(BLUP_Ind)
#########################################################
#GANHOS DE SELECAO 
##########################################################

NIS =50 #Numero de individuos selecionados


#Diferencial de selecao
Ds = mean(BLUP_Ind[1:NIS,5])-media
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
  "resultados_exemploTPDAP1planta.xlsx"
)

