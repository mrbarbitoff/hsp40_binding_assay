setwd("/media/barbitoff/DATA/Working issues/FizGen/Hsp40 mechanisms in vitro/Paper/FEMS_Submission/to_github")
library(ggplot2)
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

# Change to 'data_nm.csv'
data_q1 = read.table('data_merged.csv', sep=',', header=T)
head(data_q1)

data_q1$bound = data_q1$Pellet / (data_q1$Super + data_q1$Pellet)

ggplot(data_q1, aes(x=Ligand, y=bound, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline1 = apply(data_q1, 1, function(x) mean(data_q1[data_q1$Chaperone == x[5] & data_q1$Ligand == 0, 'bound']))
baseline1

data_q1$corr = ifelse(data_q1$Ligand == 0, 0, (data_q1$bound - baseline1)/(1-baseline1))
data_q1$corr[data_q1$corr < 0] = 0

myformula = formula(corr~Ligand/(K+Ligand))

bestfit_sis1_nm <- nls(myformula, data_q1[data_q1$Chaperone == 'Sis1', ], start=list(K=2))
summary(bestfit_sis1_nm)
bestfit_dd_nm <- nls(myformula, data_q1[data_q1$Chaperone == 'deltaDD', ], start=list(K=2))
summary(bestfit_dd_nm)
bestfit_ydj_nm <- nls(myformula, data_q1[data_q1$Chaperone == 'Ydj1', ], start=list(K=2))
summary(bestfit_ydj_nm)

coeffs1 = data.frame(K = c(summary(bestfit_sis1_nm)$coeff[1], summary(bestfit_dd_nm)$coeff[1],
                    summary(bestfit_ydj_nm)$coeff[1]), 
                    se = c(summary(bestfit_sis1_nm)$coeff[2], summary(bestfit_dd_nm)$coeff[2],
                           summary(bestfit_ydj_nm)$coeff[2]),
                    Chaperone = c('Sis1', 'deltaDD', 'Ydj1'))
coeffs1$Fibril = 'NM'

ggplot(data_q1, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol3))

ggplot(coeffs1, aes(x=Chaperone, y=K, fill=Chaperone)) + geom_bar(col='black', stat='identity') +
  geom_errorbar(ymin=coeffs1$K-coeffs1$se, ymax=coeffs1$K+coeffs1$se, width=0.5) +
  scale_y_continuous(limits=c(0, 120)) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + theme_bw() +
  ylab('Estimated Kd (uM)')


# Rnq1 goes here

# Change tp 'data_rnq.csv'
data_q2 = read.table('data_rnq1_iter2.csv', sep=',', header=T)
head(data_q2)

data_q2$bound = data_q2$Pellet / (data_q2$Super + data_q2$Pellet)

ggplot(data_q2, aes(x=Ligand, y=bound, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline2 = apply(data_q2, 1, function(x) mean(data_q2[data_q2$Chaperone == x[5] & data_q2$Ligand == 0, 'bound']))
baseline2

data_q2$corr = ifelse(data_q2$Ligand == 0, 0, (data_q2$bound - baseline2)/(1-baseline2))
data_q2$corr[data_q2$corr < 0] = 0

myformula = formula(corr~Ligand/(K+Ligand))
data_q2[data_q2$corr == 0 & data_q2$Ligand > 0, 'corr'] = 0.001


bestfit_sis1_rnq <- nls(myformula, data_q2[data_q2$Chaperone == 'Sis1', ], start=list(K=2))
summary(bestfit_sis1_rnq)
bestfit_dd_rnq <- nls(myformula, data_q2[data_q2$Chaperone == 'deltaDD', ], start=list(K=2))
summary(bestfit_dd_rnq)
#bestfit_ydj_rnq <- nls(myformula, data_q2[data_q2$Chaperone == 'Ydj1', ], start=list(K=2))
#summary(bestfit_ydj_rnq)

coeffs2 = data.frame(K = c(summary(bestfit_sis1_rnq)$coeff[1], summary(bestfit_dd_rnq)$coeff[1]), 
                    se = c(summary(bestfit_sis1_rnq)$coeff[2], summary(bestfit_dd_rnq)$coeff[2]),
                    Chaperone = c('Sis1', 'deltaDD'))
coeffs2$Fibril = 'Rnq1'

ggplot(data_q2, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol3))

ggplot(coeffs2, aes(x=Chaperone, y=K, fill=Chaperone)) + geom_bar(col='black', stat='identity') +
  geom_errorbar(ymin=coeffs2$K-coeffs2$se, ymax=coeffs2$K+coeffs2$se, width=0.5) +
  scale_y_continuous(limits=c(0, 60000)) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + theme_bw() +
  ylab('Estimated Kd (uM)')


# Hsp40 binding vs amyloidpgenic proteins
data_all= rbind(data_q1, data_q2)
data_all[data_all$corr == 0 & data_all$Ligand > 0, 'corr'] = 0.001

ggplot(data_all, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2),
                                                  control=nls.control(maxiter=100)),
                                 se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol6)) + facet_wrap(~Fibril, nrow=2)

coeffs_all = rbind(coeffs1, coeffs2)
coeffs_all

ggplot(coeffs_all, aes(x=Chaperone, y=K, fill=Chaperone)) + 
  geom_bar(col='black', stat='identity', position='dodge') +
  geom_errorbar(aes(ymin=coeffs_all$K-2*coeffs_all$se,  
                    ymax=coeffs_all$K+2*coeffs_all$se), 
                position=position_dodge(width=0.9), width=0.5) +
  scale_y_log10() +
  scale_fill_manual(values=c(mycol2, mycol1, mycol6)) + theme_bw() +
  ylab('Estimated Kd (uM)') + facet_wrap(~Fibril)

# Comparison

reg.t.test <- function(k.a, k.b, se.a, se.b, df) {
  z = (k.a - k.b)/sqrt(se.a^2 + se.b^2)
  return(pt(abs(z), lower.tail=F, df = df))
}

js = c('Sis1', 'Ydj1', 'deltaDD')
fibrils = c('NM', 'Rnq1')
# 
# for (i in fibrils) {
#   for (j in js) {
#     for (k in js) {
#       if (k != j) {
#         print(i)
#         print(j)
#         print(k)
#         print(reg.z.test(coeffs_all$K[coeffs_all$Chaperone == j & coeffs_all$Fibril == i],
#                          coeffs_all$K[coeffs_all$Chaperone == k & coeffs_all$Fibril == i],
#                          coeffs_all$se[coeffs_all$Chaperone == j & coeffs_all$Fibril == i],
#                          coeffs_all$se[coeffs_all$Chaperone == k & coeffs_all$Fibril == i]))
#       }
#     }
#   }
# }