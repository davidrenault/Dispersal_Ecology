library(GrapheR)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(dplyr)
library(vegan)
library(usdm)
library(rcompanion) 
library(FSA) 
library(ggplot2)
library(lasso2)
library('caret')
library(rcompanion)
library(car)
library(usdm)

data=read.table("data.txt", h=T, stringsAsFactors = T)
data$scale=as.factor(data$scale)
data$loc=data$LOC
data$Pgi_Q=data$Pgi/data$Q
data$Pyru_Q=data$Pyru/data$Q


###### Figure 1

p1=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=PathLength, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Mean path length (m)",x="", y = "Length (m)", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t1=p1+stat_compare_means(method = "anova") 


p2=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=S2, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Corrected sinuosity index (Benhamou, 2004)",x="", y = "S2", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t2=p2+stat_compare_means(method = "anova")



p3=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=TE, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Trajectory straightness index, E-max (Cheung et al., 2007)",x="", y = "E-max", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t3=p3+stat_compare_means(method = "anova") 



p4=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=TS, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Straightness of a trajectory (Batschelet, 1981)",x="", y = "TS", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t4=p4+stat_compare_means(method = "anova") 


p5=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=DC, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Mean directional change (DC) for all turning angles (Kitamura & Imafuku, 2015)",x="", y = "DC", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t5=p5+stat_compare_means(method = "anova") 


p6=ggplot(data, aes(x=factor(loc, levels=c("POC", "PELIS", "ISB", "VST")), y=SDDC, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot()+
  labs(title="Standard deviation of directional change (SDDC) for all turning angles (Kitamura & Imafuku, 2015)",x="", y = "SDDC", fill="locations")+
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(size=8)) 
t6=p6+stat_compare_means(method = "anova") 


ggarrange(t1, t2, t3, t4, t5, t6, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)




##### Comparison tests


############ muL
#p+stat_compare_means(method = "anova")
p+stat_compare_means(method = "kruskal.test") 

stat.test <- aov(muL ~ loc, data = data) %>% tukey_hsd()
stat.test

# general comparison between groups without parametric data
kruskal.test(muL ~ loc, data = data)

# multiple comparisons between groups
pairwise.wilcox.test(data$muL, data$loc, p.adjust.method = "BH")

dunnTest(muL ~ loc, data=data, method="bh") 


############ pathlength
#p+stat_compare_means(method = "anova")
p1+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(PathLength ~ loc, data = data)

stat.test <- aov(PathLength ~ loc, data = data) %>% tukey_hsd()
stat.test


# multiple comparisons between groups
pairwise.wilcox.test(data$PathLength, data$loc, p.adjust.method = "BH")

dunnTest(PathLength ~ loc, data=data, method="bh") 



############ S2
#p+stat_compare_means(method = "anova")
p2+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(S2 ~ loc, data = data)

stat.test <- aov(S2 ~ loc, data = data) %>% tukey_hsd()
stat.test

# multiple comparisons between groups
pairwise.wilcox.test(data$S2, data$loc, p.adjust.method = "BH")

dunnTest(S2 ~ loc, data=data, method="bh") 


############ Emax
#p3+stat_compare_means(method = "anova")
p3+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(TE ~ loc, data = data)

stat.test <- aov(TE ~ loc, data = data) %>% tukey_hsd()
stat.test


# multiple comparisons between groups
pairwise.wilcox.test(data$TE, data$loc, p.adjust.method = "BH")

dunnTest(TE ~ loc, data=data, method="bh") 


############ TS
#p4+stat_compare_means(method = "anova")
p4+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(TS ~ loc, data = data)

# multiple comparisons between groups
pairwise.wilcox.test(data$TS, data$loc, p.adjust.method = "BH")

dunnTest(TS ~ loc, data=data, method="bh") 


############ DC
#p5+stat_compare_means(method = "anova")
p5+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(DC ~ loc, data = data)

# multiple comparisons between groups
pairwise.wilcox.test(data$DC, data$loc, p.adjust.method = "BH")

dunnTest(DC ~ loc, data=data, method="bh") 

stat.test <- aov(DC ~ loc, data = data) %>% tukey_hsd()
stat.test


############ SDDC
#p6+stat_compare_means(method = "anova")
p6+stat_compare_means(method = "kruskal.test") 

# general comparison between groups without parametric data
kruskal.test(SDDC ~ loc, data = data)

# multiple comparisons between groups
pairwise.wilcox.test(data$SDDC, data$loc, p.adjust.method = "BH")

dunnTest(SDDC ~ loc, data=data, method="bh") 




###### Figure 2


# INTOC
p1=ggplot(data, aes(x=loc, y=INTOC, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "INTOC (mm)", fill="locations")+
  theme_classic()
t1=p1+stat_compare_means(method = "anova") 

# PRONO_l
p2=ggplot(data, aes(x=loc, y=PRONO_l, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "PRONO_l (mm)", fill="locations")+
  theme_classic()
t3=p2+stat_compare_means(method = "anova") 

# PRONO_w
p3=ggplot(data, aes(x=loc, y=PRONO_w, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "PRONO_w (mm)", fill="locations")+
  theme_classic()
t5=p3+stat_compare_means(method = "anova") 

# ELYT
p4=ggplot(data, aes(x=loc, y=ELYT, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "ELYT (mm)", fill="locations")+
  theme_classic()
t7=p4+stat_compare_means(method = "anova")


# STER
p5=ggplot(data, aes(x=loc, y=STER, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "STER (mm)", fill="locations")+
  theme_classic()
t9=p5+stat_compare_means(method = "anova") 

# FEMU
p6=ggplot(data, aes(x=loc, y=FEMU, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "FEMU (mm)", fill="locations")+
  theme_classic()
t11=p6+stat_compare_means(method = "anova") 


ggarrange(t1, t3, t5, t7, t9, t11, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)




###### Figure 3


# Q

p1=ggplot(data, aes(x=loc, y=Q, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "Q", fill="locations")+
  theme_classic()
t1=p1+stat_compare_means(method = "anova") # Add global p-value


# Pgi
p2=ggplot(data, aes(x=loc, y=Pgi, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "Pgi", fill="locations")+
  theme_classic()
t3=p2+stat_compare_means(method = "anova")


# Pyru
p3=ggplot(data, aes(x=loc, y=Pyru, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "Pyru", fill="locations")+
  theme_classic()
t5=p3+stat_compare_means(method = "anova")


# Pgi_Q
p4=ggplot(data, aes(x=loc, y=Pgi_Q, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "Pgi_Q", fill="locations")+
  theme_classic()
t7=p4+stat_compare_means(method = "anova")


# Pyru_Q
p5=ggplot(data, aes(x=loc, y=Pyru_Q, fill = factor(loc, levels=c("POC", "PELIS", "ISB", "VST")))) +
  geom_boxplot(show.legend = FALSE)+
  labs(x="locations", y = "Pyru_Q", fill="locations")+
  theme_classic()
t9=p5+stat_compare_means(method = "anova") 


ggarrange(t1, t3, t5, t7, t9, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)




####################
# GLM

tab <- subset(data, data$DC != 'NA') 

# Path Length
shapiro.test(tab$PathLength)
mod=glm(PathLength~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='gaussian')

step(mod)
selmod=glm(formula = PathLength ~ LOC + SEX + FEMU + STER + Pgi + Pyru + 
             LOC:SEX + SEX:STER + LOC:FEMU + LOC:STER + LOC:Pgi + LOC:Pyru, 
           family = "gaussian", data = tab)

summary(selmod)
Anova(selmod, type="III")
summary(aov(selmod))
nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)


####################
# S2
shapiro.test(tab$S2)
mod=glm(S2~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='gaussian')

step(mod)
selmod=glm(formula = S2 ~ LOC + SEX + FEMU + STER + Pgi + Pyru + LOC:SEX + 
             LOC:FEMU + LOC:STER + LOC:Pgi + LOC:Pyru, family = "gaussian", 
           data = tab)
summary(selmod)
Anova(selmod, type="III")
summary(aov(mod))
nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)




####################
# TE
shapiro.test(tab$TE)
hist(tab$TE)
mod=glm(TE~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='poisson')

step(mod)
selmod=glm(formula = TE ~ LOC * SEX + FEMU * SEX + STER * SEX + FEMU * 
             LOC + STER * LOC + LOC * Pgi + LOC * Pyru, family = "poisson", 
           data = tab)
summary(selmod)
Anova(selmod, type="III")
summary(aov(mod))
nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)

####################
# TS
shapiro.test(tab$TS)
hist(tab$TS)
mod=glm(TS~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='gaussian')

step(mod)
selmod=glm(formula = TS ~ LOC + SEX + FEMU + STER + Pgi + Pyru + SEX:FEMU + SEX:STER + 
             LOC:STER + LOC:Pyru, family = "gaussian", 
           data = tab)
summary(selmod)
Anova(selmod, type="III")
summary(aov(mod))
nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)


####################
# DC
shapiro.test(tab$DC)
hist(tab$DC)
mod=glm(DC~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='gaussian')

step(mod)
selmod=glm(formula = DC ~ LOC + SEX + FEMU + STER + Pgi + Pyru + LOC:FEMU + 
             LOC:STER + LOC:Pgi + LOC:Pyru, family = "gaussian", data = tab)
summary(selmod)
Anova(selmod, type="III")
summary(aov(mod))
nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)

####################
# SDDC
shapiro.test(tab$SDDC)
hist(tab$SDDC)
mod=glm(SDDC~LOC*SEX+FEMU*SEX+STER*SEX+FEMU*LOC+STER*LOC+LOC*Pgi+LOC*Pyru, data=tab, family='gaussian')

step(mod)
selmod=glm(formula = SDDC ~ LOC + SEX + FEMU + Pgi + Pyru + LOC:FEMU + 
             LOC:Pgi + LOC:Pyru, family = "gaussian", data = tab)
summary(selmod)
Anova(selmod, type="III")
summary(aov(mod))

nagelkerke(selmod)
library(modEvA);modEvA::Dsquared(selmod)
Dsquared(selmod, adjust = TRUE)
