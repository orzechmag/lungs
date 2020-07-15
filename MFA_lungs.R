library(FactoMineR)
library(factoextra)

notch = read.table("~/notch_mfa.txt", header = T, sep = "\t", dec = ".")
hedgehog = read.table("~/hedgehog_mfa.txt", header = T, sep = "\t", dec = ".")
erbb = read.table("~/erbb_mfa.txt", sep = "\t", dec = ".", header = T)
wnt = read.table("~/wnt_mfa.txt", sep = "\t", dec = ".", header = T)

res_notch = MFA(notch, group = c(1,1,1,1,1,1,1,174,562,144), type = c("s", "n", "n", "n", "n", "n", "n", "s", "s", "s"), ncp = 5, 
          name.group = c("age", "cancer_status", "vital_status", "stage", "smoking_history", "gender", "type", 
                         "expression_b", "expression_t", "expression_y"), num.group.sup = c(1,2,3,4,5,6,7), graph = TRUE)

res_hedgehog = MFA(hedgehog, group = c(1,1,1,1,1,1,1,210,361,37), type = c("s", "n", "n", "n", "s", "n", "n", "s", "s", "s"), ncp = 5, 
          name.group = c("age", "cancer_status", "vital_status", "stage", "smoking_history", "gender", "type", 
                         "expression_b", "expression_bl", "expression_p"), num.group.sup = c(1,2,3,4,5,6,7), graph = TRUE)


res_erbb = MFA(erbb, group = c(rep(1, 7), 86, 38, 820), type = c("s", rep("n", 6), rep("s", 3)), name.group = c("age", "relapse", "death", "stage", "smoking", 
                                                                                                                "gender", "type", "red", "purple", "turquoise"),
               num.group.sup = c(1, 2, 3, 4, 5, 6, 7, 8), graph = T)

res_wnt = MFA(wnt, group = c(rep(1, 7), 493, 98, 109), type = c("s", rep("n", 6), rep("s", 3)), name.group = c("age", "relapse", "death", "stage", "smoking", 
                                                                                                               "gender", "type", "blue", "green", "yellow"),
              num.group.sup = c(1, 2, 3, 4, 5, 6, 7, 10), graph = T)

plot.MFA(res_wnt, choix = "ind", lab.var = F, lab.ind = F, habillage = 7)

plot.MFA(res_erbb, choix = "ind", lab.var = F, lab.ind = F, habillage = 7)

plot.MFA(res_notch, axes = c(1,2), choix = "ind", new.plot = TRUE, lab.ind = FALSE, lab.par = TRUE, lab.var = TRUE, habillage = 7)

plot.MFA(res_hedgehog, axes = c(1,2), choix = "ind", new.plot = TRUE, lab.ind = FALSE, lab.par = TRUE, lab.var = TRUE, habillage = 7)
