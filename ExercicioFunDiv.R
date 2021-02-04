# Atividade diversidade funcional formigas Parque Nacional Catimbau - Gabriel Boldorini ####

# Pacotes ####
install.packages("FD")
install.packages("dplyr")
install.packages("PerformanceAnalytics")
install.packages("ggplot2")
library("FD")
library("dplyr")
library("PerformanceAnalytics")
library("ggplot2")

# Diret?rio ####
setwd("C:/Users/joyce/OneDrive/Documentos/Gabriel/Diversidade funcional/Atividade R - FunDiv")

# Dados ####
sp_traits <- read.table("tra?os_esp?cies_catimbau.txt", h=T, row.names = 1, sep="\t")
sp_traits <- sp_traits[,1:18] # remo??o de coluna que foi inserida na entrada dos dados no r
comunidades <- read.table("comunidades_formigas_catimbau.txt", h=T, row.names = 1, sep="\t")
metadata_catimbau <-  read.table("metadata_catimbau.txt", h=T, row.names = 1, sep="\t")

# Traits ####
traits_hl_p_cv <- select(sp_traits, "xHL","xP", "xCV")

# Verificar correla??o entre os traits selecionados ####
chart.Correlation(traits_hl_p_cv, histogram=TRUE, pch=19)

# Verificar a distribui??o dos traits ####
par(mfrow = c(1,3), cex = 0.5)
head_length <- density(traits_hl_p_cv$xHL, na.rm = T)
plot(head_length, main = "Head length")
pilosity <- density(traits_hl_p_cv$xP, na.rm = T)
plot(pilosity, main = "Pilosity")
colour_value <- density(traits_hl_p_cv$xCV, na.rm = T)
plot(colour_value, main = "Colour Value")

# Community Weigthed Mean ####
functcomp(x=traits_hl_p_cv, a=as.matrix(t(comunidades)), CWM.type = c("dom"))
functcomp(x=traits_hl_p_cv, a=as.matrix(t(comunidades)), CWM.type = c("all"))

# dbFD ####
rownames(traits_hl_p_cv)==rownames(comunidades)
div.func <- dbFD(x=traits_hl_p_cv, a=t(comunidades), corr = "lingoes")

# Ind?ces de diversidade funcional ####
FRic <- div.func$FRic
Feve <- div.func$FEve
FDiv <- div.func$FDiv
FDis <- div.func$FDis
RaoQ <- div.func$RaoQ

# CWM para cada trait ####
cwm_head <- div.func$CWM$xHL
cwm_pilosity <- div.func$CWM$xP
cwm_colour <- div.func$CWM$xCV

# Juntar os dados em uma ?nica tabela ####
dados_funcionais <- cbind(FRic,Feve,FDiv,FDis,RaoQ, cwm_head, cwm_pilosity, cwm_colour)
dados_completos <- cbind(metadata_catimbau, dados_funcionais)

# Ind?ces de diversidade para a comunidade ####
abundancia_comunidade <- apply(t(comunidades), 1, sum)
H <- diversity (t(comunidades), index = "shannon")
S <- specnumber (t(comunidades))
J <- H/log(specnumber(t(comunidades)))

# Tabela final ####
dados_completos2<-cbind(dados_completos,abundancia_comunidade,H,S,J)

# Modelos - Ind?ces de diversidade funcional ####
model_FRic <- lm(FRic ~ GMDI + DH, data=dados_completos2)
summary(model_FRic) # valores significativos e negativo para perturba??o e positivo para d?ficit h?drico

model_Feve <- lm(Feve ~ GMDI + DH, data=dados_completos2)
summary(model_Feve)

model_RaoQ <- lm(RaoQ ~ GMDI + DH, data=dados_completos2)
summary(model_RaoQ) # valor significativo e negativo para d?ficit h?drico

# Plotar o gr?fico de riqueza funcional e d?fict h?drico e perturba??o ####
ggplot(dados_completos2,aes(DH, FRic)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="D?ficit H?drico", y="Riqueza funcional")

ggplot(dados_completos2,aes(GMDI, FRic)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Perturba??o", y="Riqueza funcional")

# Plotar o gr?fico de dissimilaridade funcional e d?fict h?drico ####
ggplot(dados_completos2,aes(DH, RaoQ)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="D?ficit H?drico", y="Dissimilaridade funcional")


# Modelos - CWM para cada trait ####
model_head <- lm(cwm_head ~ GMDI + DH, data=dados_completos2)
summary(model_head) # valor significativo e negativo para d?ficit h?drico

model_pilosity <-lm(cwm_pilosity ~ GMDI + DH, data=dados_completos2)
summary(model_pilosity) # valor significativo e negativo para perturba??o

model_colour <- lm(cwm_colour ~ GMDI + DH, data=dados_completos2)
summary(model_colour) # valor significativo e positivo para d?ficit h?drico

# Plotar o gr?fico do CWM do trait Head length e d?ficit h?drico ####
ggplot(dados_completos2,aes(DH, cwm_head)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="D?ficit h?drico", y="Comprimento da cabe?a")

# Plotar o gr?fico do CWM do trait Pilosity e perturba??o ####
ggplot(dados_completos2,aes(GMDI, cwm_pilosity)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Perturba??o", y="Pilosidade")

# Plotar o gr?fico do CWM do trait Colour value e d?ficit h?drico ####
ggplot(dados_completos2,aes(DH, cwm_colour)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="D?ficit h?drico", y="Valor de colora??o")


# Respostas da atividade ####
# 1)	Indique e selecione os tra?os morfol?gicos que voc? acredita que podem ser importantes para explicar as
# respostas das formigas ?s mudan?as ambientais, especialmente ?s mudan?as na intensidade da perturba??o e de
# aridez. Justifique sua resposta.

# R: Selecionei os traits de comprimento da cabe?a, pilosidade e colora??o. O primeiro foi
# ecolhido por indicar a dieta do animal, pois esp?cies com diferentes dietas podem ser afetadas de diferentes
# formas pelas perturba??es e aridez (mudan?a de recursos). O segundo trait foi escolhido por refletir uma
# toler?ncia a desidrata??o, portanto em diferentes gradientes de aridez isso pode ser importante. O terceiro
# trait foi selecionado pois pode indicar tanto uma resposta termal quanto a outros estresses ambientais.

# 2)	Quais dos tra?os selecionados respondem aos gradientes de perturba??o e aridez? Os tra?os que respondem, em
# que sentido ? a resposta ao longo desses gradientes? Explique detalhadamente todos os passos que seguiu para
# chegar ? sua resposta.

# R: Calculei diferentes ?ndices de diversidade funcional para cada comunidade (com os 3 traits selecionados)
# e o cwm para cada trait em cada comunidade. Os valores significativos foram para os modelos de riqueza e
# dissimilaridade funcional. Onde para o modelo FRic foi encontrado uma rela??o negativa para perturba??o e
# positiva para aridez. No modelo de RaoQ foi encontrado valor significativo e negativo apenas para aridez.
# Al?m disso, para os modelos de cwm de cada trait, foi encontrado rela??o significativa e negativa entre o
# comprimento da cabe?a e o d?ficit h?drico; significativa e negativa entre pilosidade e perturba??o;
# significativa e positiva entre a colora??o e o d?ficit h?drico.
